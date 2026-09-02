#!/usr/bin/env python
"""Plot omega_R along the adiabatic spin-down trajectory, with backreaction
on (evolving Leaver frequency, M_BH and chi both respond to the extracted
cloud) versus off (frequency frozen at its initial value, so the trajectory
is the exact-solvable straight-line J_BH(M_BH) control already validated in
adiabatic_endpoint.validate_fixed_frequency_trajectory).

Both trajectories start at the same (M_BH, chi, alpha) and integrate the
horizon-flux ODE dJ_BH/dM_BH = m/omega_R to the superradiant threshold
omega_R = m*Omega_H(M_BH, chi). "Backreaction off" freezes omega_R in that
ODE at its initial value (the cloud still drains the hole, but the frequency
never updates -- the naïve control used elsewhere in this repo).
"Backreaction on" re-solves the Leaver problem at every step, so the
spectral feedback from the shrinking mass/spin is included.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from adiabatic_endpoint import horizon_frequency, physical_mode_frequency


def integrate_trajectory(alpha, initial_dimensionless_spin,
                         initial_black_hole_mass, ell, azimuthal_number,
                         principal_quantum_number, evolve_frequency):
    boson_mass = alpha / initial_black_hole_mass
    initial_black_hole_angular_momentum = (
        initial_dimensionless_spin * initial_black_hole_mass**2
    )
    initial_frequency = physical_mode_frequency(
        initial_black_hole_mass, initial_black_hole_angular_momentum,
        boson_mass, ell, azimuthal_number, principal_quantum_number
    ).real

    def frequency_at(black_hole_mass, black_hole_angular_momentum):
        if not evolve_frequency:
            return initial_frequency
        return physical_mode_frequency(
            black_hole_mass, black_hole_angular_momentum, boson_mass,
            ell, azimuthal_number, principal_quantum_number
        ).real

    def angular_momentum_derivative(black_hole_mass, angular_momentum_state):
        frequency = frequency_at(black_hole_mass, angular_momentum_state[0])
        return np.array((azimuthal_number / frequency,))

    def superradiant_threshold(black_hole_mass, angular_momentum_state):
        black_hole_angular_momentum = angular_momentum_state[0]
        dimensionless_spin = (
            black_hole_angular_momentum / black_hole_mass**2
        )
        frequency = frequency_at(black_hole_mass, black_hole_angular_momentum)
        return (
            frequency
            - azimuthal_number
            * horizon_frequency(black_hole_mass, dimensionless_spin)
        )

    superradiant_threshold.terminal = True
    superradiant_threshold.direction = 1.0
    solution = solve_ivp(
        angular_momentum_derivative,
        (initial_black_hole_mass, 0.5 * initial_black_hole_mass),
        np.array((initial_black_hole_angular_momentum,)),
        events=superradiant_threshold, dense_output=True,
        rtol=1.0e-10, atol=1.0e-12, max_step=0.005 * initial_black_hole_mass,
    )
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError("trajectory integration failed to reach threshold")

    final_black_hole_mass = float(solution.t_events[0][0])
    black_hole_mass_grid = np.linspace(
        initial_black_hole_mass, final_black_hole_mass, 400
    )
    black_hole_angular_momentum_grid = solution.sol(black_hole_mass_grid)[0]
    dimensionless_spin_grid = (
        black_hole_angular_momentum_grid / black_hole_mass_grid**2
    )
    frequency_grid = np.array([
        frequency_at(mass, angular_momentum) for mass, angular_momentum
        in zip(black_hole_mass_grid, black_hole_angular_momentum_grid)
    ])
    horizon_frequency_grid = horizon_frequency(
        black_hole_mass_grid, dimensionless_spin_grid
    )
    cloud_mass_fraction_grid = (
        1.0 - black_hole_mass_grid / initial_black_hole_mass
    )
    return dict(
        black_hole_mass=black_hole_mass_grid,
        dimensionless_spin=dimensionless_spin_grid,
        cloud_mass_fraction=cloud_mass_fraction_grid,
        omega_R=frequency_grid,
        omega_horizon=horizon_frequency_grid,
        final_black_hole_mass=final_black_hole_mass,
    )


def main():
    alpha = 0.42
    initial_dimensionless_spin = 0.99
    ell, azimuthal_number, principal_quantum_number = 1, 1, 2

    trajectory_off = integrate_trajectory(
        alpha, initial_dimensionless_spin, 1.0, ell, azimuthal_number,
        principal_quantum_number, evolve_frequency=False
    )
    trajectory_on = integrate_trajectory(
        alpha, initial_dimensionless_spin, 1.0, ell, azimuthal_number,
        principal_quantum_number, evolve_frequency=True
    )
    print(f"backreaction OFF: final M_BH = {trajectory_off['final_black_hole_mass']:.8f}")
    print(f"backreaction ON:  final M_BH = {trajectory_on['final_black_hole_mass']:.8f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    ax[0].plot(
        trajectory_off["cloud_mass_fraction"], trajectory_off["omega_R"],
        color="C1", lw=2, label=r"backreaction OFF ($\omega_R$ frozen)"
    )
    ax[0].plot(
        trajectory_on["cloud_mass_fraction"], trajectory_on["omega_R"],
        color="C0", lw=2, label=r"backreaction ON (Leaver $\omega_R$ re-solved)"
    )
    ax[0].plot(
        trajectory_on["cloud_mass_fraction"],
        azimuthal_number * trajectory_on["omega_horizon"],
        color="k", ls="--", lw=1, label=r"$m\,\Omega_H(t)$ (moving threshold)"
    )
    ax[0].axhline(
        azimuthal_number * trajectory_off["omega_horizon"][0], color="0.6",
        ls=":", lw=1, label=r"$m\,\Omega_H(0)$ (fixed threshold)"
    )
    ax[0].set_xlabel(r"cloud mass fraction, $1 - M_{\rm BH}/M_{\rm BH,0}$")
    ax[0].set_ylabel(r"$\omega_R$")
    ax[0].set_title(rf"$\alpha={alpha}$, $\chi_0={initial_dimensionless_spin}$, mode $(n,\ell,m)=({principal_quantum_number},{ell},{azimuthal_number})$")
    ax[0].legend(fontsize=8)

    ax[1].plot(
        trajectory_off["cloud_mass_fraction"], trajectory_off["dimensionless_spin"],
        color="C1", lw=2, label="backreaction OFF"
    )
    ax[1].plot(
        trajectory_on["cloud_mass_fraction"], trajectory_on["dimensionless_spin"],
        color="C0", lw=2, label="backreaction ON"
    )
    ax[1].set_xlabel(r"cloud mass fraction, $1 - M_{\rm BH}/M_{\rm BH,0}$")
    ax[1].set_ylabel(r"Kerr spin $\chi = J_{\rm BH}/M_{\rm BH}^2$")
    ax[1].set_title("spin-down trajectory to threshold")
    ax[1].legend(fontsize=8)

    fig.tight_layout()
    output_path = os.path.join(_HERE, "omega_trajectory.png")
    fig.savefig(output_path, dpi=130)
    print("saved", output_path)


if __name__ == "__main__":
    main()
