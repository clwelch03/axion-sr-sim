"""Time-domain growth-rate measurement following Dolan 1212.1477.

This module will evolve broadband Gaussian initial data and extract individual
quasi-bound-state growth rates with the frequency filter in
``frequency_filter.py``.  It replaces the earlier eigenmode-seed/log-norm
driver, whose fitted rate was PML-dependent and was not Dolan's Fig. 13 method.
"""

import argparse
import os
import time

import numpy as np

from frequency_filter import (
    find_mode_frequency,
    fit_growth_rate,
    fourier_amplitudes,
    reconstruct_filtered_mode,
)


def dolan_gaussian_initial_state(tower, gaussian_center=20.0,
                                 gaussian_width=10.0):
    """Return Dolan's time-symmetric massive-field Gaussian initial state.

    Equation (45) defines the Gaussian with psi_dot=0 and only the lowest
    angular channel populated.  Figure 11 uses x0=20M and sigma_x=10M for the
    massive-field example; using those values for Fig. 13 is a documented
    inference because its caption does not repeat the initial-data parameters.
    """
    if not np.isfinite(gaussian_center):
        raise ValueError("gaussian_center must be finite")
    if not np.isfinite(gaussian_width) or gaussian_width <= 0.0:
        raise ValueError("gaussian_width must be finite and positive")

    psi = np.zeros(
        (tower.n_channels, tower.num_points), dtype=complex)
    psi[0] = np.exp(
        -0.5 * ((tower.x - gaussian_center) / gaussian_width) ** 2)
    pi = np.zeros_like(psi)
    chi = np.zeros_like(psi)

    # Tower.rhs freezes two nodes at each hard boundary; impose the Dirichlet
    # data exactly in the initial state as well.
    psi[:, :2] = 0.0
    psi[:, -2:] = 0.0
    return psi, pi, chi


def evolve_tower_and_sample_lowest_channel(
        tower, initial_state, time_final, probe_x, sample_every_steps=1,
        courant_factor=0.8, checkpoint_path=None,
        checkpoint_every_steps=None, resume=False):
    """Evolve and uniformly sample the lowest-l field at one fixed x.

    Dolan uses ``dt = 0.8 dx`` and Fourier transforms the complex time series
    of each angular channel at a fixed tortoise coordinate.  The requested
    probe position is snapped to the nearest grid point; its actual coordinate
    is returned so that this discretization is visible to the caller.
    """
    if not np.isfinite(time_final) or time_final <= 0.0:
        raise ValueError("time_final must be finite and positive")
    if not np.isfinite(probe_x):
        raise ValueError("probe_x must be finite")
    if (not isinstance(sample_every_steps, (int, np.integer))
            or sample_every_steps <= 0):
        raise ValueError("sample_every_steps must be a positive integer")
    if not np.isfinite(courant_factor) or courant_factor <= 0.0:
        raise ValueError("courant_factor must be finite and positive")

    time_step = courant_factor * tower.dx
    total_steps = round(time_final / time_step)
    if not np.isclose(total_steps * time_step, time_final,
                      rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("time_final must be an integer multiple of dt")
    if total_steps % sample_every_steps != 0:
        raise ValueError(
            "total time steps must be divisible by sample_every_steps")

    if checkpoint_path is not None:
        if (not isinstance(checkpoint_every_steps, (int, np.integer))
                or checkpoint_every_steps <= 0):
            raise ValueError(
                "checkpoint_every_steps must be a positive integer")
        if checkpoint_every_steps % sample_every_steps != 0:
            raise ValueError(
                "checkpoint_every_steps must be divisible by sample_every_steps")
    elif resume:
        raise ValueError("resume=True requires checkpoint_path")

    probe_index = int(np.argmin(np.abs(tower.x - probe_x)))
    actual_probe_x = float(tower.x[probe_index])
    number_of_samples = total_steps // sample_every_steps + 1
    sample_times = np.arange(number_of_samples, dtype=float)
    sample_times *= sample_every_steps * time_step
    probe_signal = np.empty(number_of_samples, dtype=complex)

    if resume:
        with np.load(checkpoint_path) as checkpoint:
            checkpoint_parameters = np.array([
                tower.spin,
                tower.alpha,
                tower.m,
                tower.n,
                tower.x_pml_center,
                tower.pml_width,
            ])
            if not np.array_equal(
                    checkpoint["tower_parameters"], checkpoint_parameters):
                raise ValueError("checkpoint parameters do not match tower")
            if not np.array_equal(checkpoint["grid_x"], tower.x):
                raise ValueError("checkpoint grid does not match tower")
            if not np.array_equal(
                    checkpoint["l_values"], np.asarray(tower.l_values)):
                raise ValueError("checkpoint channels do not match tower")
            if not np.isclose(checkpoint["time_step"], time_step):
                raise ValueError("checkpoint time step does not match")
            if int(checkpoint["sample_every_steps"]) != sample_every_steps:
                raise ValueError("checkpoint sampling cadence does not match")
            if not np.isclose(checkpoint["actual_probe_x"], actual_probe_x):
                raise ValueError("checkpoint probe position does not match")

            completed_steps = int(checkpoint["completed_steps"])
            if completed_steps > total_steps:
                raise ValueError("checkpoint is later than requested time_final")
            state = tuple(
                np.array(checkpoint[field_name], copy=True)
                for field_name in ("psi", "pi", "chi")
            )
            saved_probe_signal = checkpoint["probe_signal"]
            saved_sample_times = checkpoint["sample_times"]
            sample_index = saved_probe_signal.size
            expected_sample_index = completed_steps // sample_every_steps + 1
            if sample_index != expected_sample_index:
                raise ValueError("checkpoint sample count is inconsistent")
            if not np.array_equal(
                    saved_sample_times, sample_times[:sample_index]):
                raise ValueError("checkpoint sample times are inconsistent")
            probe_signal[:sample_index] = saved_probe_signal
    else:
        completed_steps = 0
        state = tuple(np.array(field, copy=True) for field in initial_state)
        probe_signal[0] = state[0][0, probe_index]
        sample_index = 1

    def save_checkpoint(step_index):
        temporary_path = os.fspath(checkpoint_path) + ".tmp"
        with open(temporary_path, "wb") as checkpoint_file:
            np.savez(
                checkpoint_file,
                tower_parameters=np.array([
                    tower.spin,
                    tower.alpha,
                    tower.m,
                    tower.n,
                    tower.x_pml_center,
                    tower.pml_width,
                ]),
                grid_x=tower.x,
                l_values=np.asarray(tower.l_values),
                time_step=time_step,
                sample_every_steps=sample_every_steps,
                actual_probe_x=actual_probe_x,
                completed_steps=step_index,
                sample_times=sample_times[:sample_index],
                probe_signal=probe_signal[:sample_index],
                psi=state[0],
                pi=state[1],
                chi=state[2],
            )
        os.replace(temporary_path, checkpoint_path)

    for step_index in range(completed_steps + 1, total_steps + 1):
        state = tower.rk4_step(state, time_step)
        if step_index % sample_every_steps == 0:
            probe_signal[sample_index] = state[0][0, probe_index]
            sample_index += 1
        if (checkpoint_path is not None
                and step_index % checkpoint_every_steps == 0):
            save_checkpoint(step_index)

    if not all(np.all(np.isfinite(field)) for field in state):
        raise FloatingPointError("evolution produced a non-finite state")
    if checkpoint_path is not None and total_steps % checkpoint_every_steps:
        save_checkpoint(total_steps)
    return sample_times, probe_signal, actual_probe_x


def measure_filtered_mode(times, probe_signal, approximate_mode_frequency,
                          search_half_width, filter_bin_count=5.0,
                          central_fraction=0.5):
    """Locate, reconstruct, and measure one mode by Dolan's procedure.

    ``approximate_mode_frequency`` and ``search_half_width`` deliberately make
    the expected spectral neighborhood explicit.  They should come from an
    independent frequency-domain estimate, rather than an unconstrained peak
    search that could label a broadband transient as a quasi-bound state.

    Returns the interpolated physical frequency, fitted growth rate, and the
    reconstructed complex mode.  Keeping the reconstruction lets the caller
    inspect its central log-amplitude before accepting the fitted rate.
    """
    angular_frequencies, amplitudes = fourier_amplitudes(times, probe_signal)
    measured_mode_frequency = find_mode_frequency(
        angular_frequencies,
        amplitudes,
        approximate_mode_frequency,
        search_half_width,
    )
    filtered_mode = reconstruct_filtered_mode(
        times,
        probe_signal,
        measured_mode_frequency,
        filter_bin_count=filter_bin_count,
    )
    measured_growth_rate = fit_growth_rate(
        times,
        filtered_mode,
        central_fraction=central_fraction,
    )
    return measured_mode_frequency, measured_growth_rate, filtered_mode


def main():
    """Run or resume one checkpointed Dolan time-domain record."""
    parser = argparse.ArgumentParser(
        description="Evolve Dolan Gaussian data and save a fixed-x time series")
    parser.add_argument("--output", required=True,
                        help="checkpoint/result .npz path")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--time-final", type=float, default=100000.0)
    parser.add_argument("--spin", type=float, default=0.99)
    parser.add_argument("--alpha", type=float, default=0.42)
    parser.add_argument("--m", type=int, default=1)
    parser.add_argument("--l-max", type=int, default=5)
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--x-left", type=float, default=-200.0)
    parser.add_argument("--x-right", type=float, default=511.0)
    parser.add_argument("--x-pml-center", type=float, default=-100.0)
    parser.add_argument("--pml-width", type=float, default=10.0)
    parser.add_argument("--probe-x", type=float, default=20.0)
    parser.add_argument("--sample-every-steps", type=int, default=10)
    parser.add_argument("--checkpoint-interval", type=float, default=5000.0,
                        help="physical time between checkpoints")
    arguments = parser.parse_args()

    from solver import Tower

    tower = Tower(
        spin=arguments.spin,
        alpha=arguments.alpha,
        m=arguments.m,
        l_max=arguments.l_max,
        x_left=arguments.x_left,
        x_right=arguments.x_right,
        n=arguments.n,
        x_pml_center=arguments.x_pml_center,
        pml_width=arguments.pml_width,
    )
    initial_state = dolan_gaussian_initial_state(tower)
    time_step = 0.8 * tower.dx
    checkpoint_every_steps = round(
        arguments.checkpoint_interval / time_step)
    if not np.isclose(
            checkpoint_every_steps * time_step,
            arguments.checkpoint_interval,
            rtol=1.0e-12,
            atol=1.0e-12):
        raise ValueError("checkpoint_interval must be an integer multiple of dt")

    start_wall_time = time.monotonic()
    sample_times, probe_signal, actual_probe_x = (
        evolve_tower_and_sample_lowest_channel(
            tower,
            initial_state,
            time_final=arguments.time_final,
            probe_x=arguments.probe_x,
            sample_every_steps=arguments.sample_every_steps,
            checkpoint_path=arguments.output,
            checkpoint_every_steps=checkpoint_every_steps,
            resume=arguments.resume,
        )
    )
    elapsed_wall_time = time.monotonic() - start_wall_time
    print(
        f"completed t={sample_times[-1]:.1f} with {sample_times.size} samples "
        f"at x={actual_probe_x:.6f} in {elapsed_wall_time:.1f} s; "
        f"max |psi|={np.max(np.abs(probe_signal)):.6e}",
        flush=True,
    )


if __name__ == "__main__":
    main()
