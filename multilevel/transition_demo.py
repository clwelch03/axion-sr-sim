#!/usr/bin/env python
"""Cheap coherent transition demonstration for the multilevel scaffold.

This is a projected amplitude model, not a calibrated radial evolution:
horizon fluxes and backreaction are intentionally omitted.  It does include
the reduced-equation radial overlap from the Leaver profiles.
It integrates

    dA_j/dt = -i [omega_j A_j + lambda sum Q[j,a,b,c] A_a A_b* A_c]

using real Leaver frequencies so the test remains Hamiltonian.  The active
set contains |211>, |100>, and |322>; the latter is allowed through the
near-resonant combination |211> + |100> -> |322>.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from couplings import cubic_force, radial_weighted_coupling_tensor
from mode_catalog import MODE_100, MODE_211, MODE_322, leaver_frequency
from radial1d import build_grid


def run(alpha=0.30, dimensionless_spin=0.99, interaction_strength=1000.0,
        final_time=300.0, time_step=0.002, recorded_points=601,
        xmin=-40.0, xmax=120.0, grid_points=600, series_terms=400,
        initial_amplitudes=(0.8 + 0.0j, 0.03 + 0.0j, 0.0 + 0.0j),
        loss_rates=None, use_leaver_linear_rates=False):
    """Return time series for the reduced three-mode transition test.

    ``loss_rates`` is an optional real array, one entry per active mode, that
    adds an explicit non-Hermitian horizon/continuum sink ``-loss_rates * A``
    to the amplitude equation.  A positive entry depletes that channel's
    occupation; it is independent of, and not inferred from, the Hermitian
    angular/radial coupling tensor.  ``None`` (the default) disables every
    sink and keeps the bound-only system exactly Hermitian.

    ``use_leaver_linear_rates`` is an opt-in flag (default ``False``, so the
    existing exact-norm-conservation behavior of this function is unchanged
    unless a caller explicitly asks for it).  When ``True`` each channel's own
    ``Im(complex_frequency)`` from the Leaver solve is added directly as a
    physically calibrated linear growth/decay term, ``+Im(omega_j) * A_j``,
    on top of the existing cubic term.  This is not a hand-tuned rate: it is
    exactly the imaginary part of that mode's own quasi-bound-state frequency,
    positive (superradiant growth) or negative (net absorption/decay)
    whichever the Leaver solve gives.  It composes additively with any
    user-supplied ``loss_rates`` -- both sinks/sources apply simultaneously if
    both are given.
    """
    x, radius = build_grid(
        dimensionless_spin, xmin, xmax, grid_points
    )
    active_modes, radial_modes, coupling_tensor = radial_weighted_coupling_tensor(
        (MODE_211, MODE_100, MODE_322), alpha, dimensionless_spin,
        x, radius, series_terms=series_terms
    )
    complex_frequencies = np.array([
        leaver_frequency(mode_spec, alpha, dimensionless_spin)
        for mode_spec in active_modes
    ])
    frequencies = complex_frequencies.real
    amplitudes = np.asarray(initial_amplitudes, dtype=complex).copy()
    if amplitudes.shape != (len(active_modes),):
        raise ValueError("initial_amplitudes must have one entry per active mode")
    if loss_rates is None:
        loss_rates = np.zeros(len(active_modes))
    else:
        loss_rates = np.asarray(loss_rates, dtype=float)
        if loss_rates.shape != (len(active_modes),):
            raise ValueError("loss_rates must have one entry per active mode")
    if use_leaver_linear_rates:
        leaver_linear_rates = complex_frequencies.imag
    else:
        leaver_linear_rates = np.zeros(len(active_modes))

    def right_hand_side(current_amplitudes):
        return (
            -1j * (
                frequencies * current_amplitudes
                + interaction_strength
                * cubic_force(current_amplitudes, coupling_tensor)
            )
            - loss_rates * current_amplitudes
            + leaver_linear_rates * current_amplitudes
        )

    step_count = int(round(final_time / time_step))
    recording_indices = np.linspace(
        0, step_count, recorded_points, dtype=int
    )
    recorded_times = []
    recorded_amplitudes = []
    next_recording = 0
    for step_index in range(step_count + 1):
        if (next_recording < recorded_points
                and step_index >= recording_indices[next_recording]):
            recorded_times.append(step_index * time_step)
            recorded_amplitudes.append(amplitudes.copy())
            next_recording += 1
        if step_index == step_count:
            break
        first_stage = right_hand_side(amplitudes)
        second_stage = right_hand_side(
            amplitudes + 0.5 * time_step * first_stage
        )
        third_stage = right_hand_side(
            amplitudes + 0.5 * time_step * second_stage
        )
        fourth_stage = right_hand_side(
            amplitudes + time_step * third_stage
        )
        amplitudes += (time_step / 6.0) * (
            first_stage + 2.0 * second_stage
            + 2.0 * third_stage + fourth_stage
        )

    recorded_amplitudes = np.asarray(recorded_amplitudes)
    populations = np.abs(recorded_amplitudes) ** 2
    return {
        "active_modes": active_modes,
        "complex_frequencies": complex_frequencies,
        "radial_modes": radial_modes,
        "x": x,
        "radius": radius,
        "t": np.asarray(recorded_times),
        "amplitudes": recorded_amplitudes,
        "populations": populations,
        "norm": np.sum(populations, axis=1),
        "interaction_strength": interaction_strength,
    }


if __name__ == "__main__":
    output = run()
    labels = [mode.label for mode in output["active_modes"]]
    populations = output["populations"]
    print("active modes:", labels)
    print("real frequencies:", output["complex_frequencies"].real)
    print("imaginary Leaver rates (omitted in demo):",
          output["complex_frequencies"].imag)
    print("max populations:", dict(zip(labels, populations.max(axis=0))))
    print("final populations:", dict(zip(labels, populations[-1])))
    print("maximum norm drift:",
          np.max(np.abs(output["norm"] - output["norm"][0])))
