#!/usr/bin/env python
"""Backreaction-only saturation endpoint for one scalar Kerr bound level.

This is the lossless, single-level Tier-1 baseline.  The cloud is assumed to
remain in one instantaneous Leaver eigenmode while its scaled occupation N
grows adiabatically.  Each incremental horizon transfer obeys

    dM_BH / dN = -omega_R,
    dJ_BH / dN = -m,

or equivalently dJ_BH/dM_BH = m/omega_R.  The endpoint is found by integrating
this flux trajectory to the superradiant threshold omega_R = m Omega_H.  Cloud
energy and angular momentum are the cumulative losses from the black hole; the
code does not incorrectly revalue the whole cloud as omega_R changes.

No time evolution, nonlinear detuning, level transitions, gravitational
radiation, or axion radiation is included here.
"""

import argparse
from dataclasses import dataclass
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


_HERE = os.path.dirname(os.path.abspath(__file__))
_LINEAR_RATES = os.path.abspath(
    os.path.join(_HERE, "..", "..", "linear_rates")
)
sys.path.insert(0, _LINEAR_RATES)

import leaver


@dataclass(frozen=True)
class BackreactionEndpoint:
    initial_black_hole_mass: float
    initial_dimensionless_spin: float
    initial_cloud_mass_fraction: float
    initial_scaled_occupation: float
    final_black_hole_mass: float
    final_dimensionless_spin: float
    final_cloud_mass_fraction: float
    final_scaled_occupation: float
    final_frequency: float
    final_horizon_frequency: float
    fractional_horizon_area_change: float
    energy_conservation_error: float
    angular_momentum_conservation_error: float
    threshold_error: float


@dataclass(frozen=True)
class FixedFrequencyValidation:
    fixed_frequency: float
    numerical_final_black_hole_mass: float
    analytic_final_black_hole_mass: float
    black_hole_mass_error: float
    angular_momentum_error: float
    numerical_threshold_error: float


def horizon_frequency(black_hole_mass, dimensionless_spin):
    """Return Omega_H for a Kerr black hole in G=c=1 units."""
    horizon_factor = 1.0 + np.sqrt(1.0 - dimensionless_spin**2)
    return dimensionless_spin / (2.0 * black_hole_mass * horizon_factor)


def horizon_area(black_hole_mass, dimensionless_spin):
    """Return the Kerr horizon area in G=c=1 units."""
    horizon_factor = 1.0 + np.sqrt(1.0 - dimensionless_spin**2)
    return 8.0 * np.pi * black_hole_mass**2 * horizon_factor


def physical_mode_frequency(black_hole_mass, black_hole_angular_momentum,
                            boson_mass, ell, azimuthal_number,
                            principal_quantum_number):
    """Return physical omega from Leaver's dimensionless M_BH*omega."""
    dimensionless_spin = (
        black_hole_angular_momentum / black_hole_mass**2
    )
    alpha = black_hole_mass * boson_mass
    dimensionless_frequency = complex(leaver.solve(
        alpha, dimensionless_spin, l=ell, m=azimuthal_number,
        n=principal_quantum_number
    ))
    return dimensionless_frequency / black_hole_mass


def validate_fixed_frequency_trajectory(
        alpha=0.42, initial_dimensionless_spin=0.99,
        initial_black_hole_mass=1.0, ell=1, azimuthal_number=1,
        principal_quantum_number=2):
    """Compare the numerical trajectory to its constant-frequency solution.

    With omega_R held fixed, the horizon-flux equation has the exact solution

        J_BH(M_BH) = J_BH,0 + (m / omega_R) (M_BH - M_BH,0).

    Only the threshold mass requires a one-dimensional root solve.  This is an
    exact-answer test of the ODE integration and terminal-event machinery.
    """
    initial_black_hole_angular_momentum = (
        initial_dimensionless_spin * initial_black_hole_mass**2
    )
    boson_mass = alpha / initial_black_hole_mass
    fixed_frequency = physical_mode_frequency(
        initial_black_hole_mass, initial_black_hole_angular_momentum,
        boson_mass, ell, azimuthal_number, principal_quantum_number
    ).real

    def analytic_angular_momentum(black_hole_mass):
        return (
            initial_black_hole_angular_momentum
            + (azimuthal_number / fixed_frequency)
            * (black_hole_mass - initial_black_hole_mass)
        )

    def threshold_from_analytic_trajectory(black_hole_mass):
        black_hole_angular_momentum = analytic_angular_momentum(
            black_hole_mass
        )
        dimensionless_spin = (
            black_hole_angular_momentum / black_hole_mass**2
        )
        if not 0.0 < dimensionless_spin < 1.0:
            return np.nan
        return (
            fixed_frequency
            - azimuthal_number
            * horizon_frequency(black_hole_mass, dimensionless_spin)
        )

    trial_masses = np.linspace(
        initial_black_hole_mass, 0.5 * initial_black_hole_mass, 501
    )
    analytic_bracket = None
    previous_mass = trial_masses[0]
    previous_residual = threshold_from_analytic_trajectory(previous_mass)
    for trial_mass in trial_masses[1:]:
        trial_residual = threshold_from_analytic_trajectory(trial_mass)
        if np.isfinite(previous_residual) and np.isfinite(trial_residual):
            if previous_residual * trial_residual <= 0.0:
                analytic_bracket = (trial_mass, previous_mass)
                break
        previous_mass = trial_mass
        previous_residual = trial_residual
    if analytic_bracket is None:
        raise RuntimeError("could not bracket the fixed-frequency threshold")

    analytic_final_black_hole_mass = brentq(
        threshold_from_analytic_trajectory, *analytic_bracket,
        xtol=1.0e-14, rtol=1.0e-14
    )
    analytic_final_black_hole_angular_momentum = (
        analytic_angular_momentum(analytic_final_black_hole_mass)
    )

    def constant_angular_momentum_derivative(_black_hole_mass,
                                             _angular_momentum_state):
        return np.array((azimuthal_number / fixed_frequency,))

    def numerical_threshold(black_hole_mass, angular_momentum_state):
        black_hole_angular_momentum = angular_momentum_state[0]
        dimensionless_spin = (
            black_hole_angular_momentum / black_hole_mass**2
        )
        return (
            fixed_frequency
            - azimuthal_number
            * horizon_frequency(black_hole_mass, dimensionless_spin)
        )

    numerical_threshold.terminal = True
    numerical_threshold.direction = 1.0
    numerical_solution = solve_ivp(
        constant_angular_momentum_derivative,
        (initial_black_hole_mass, 0.5 * initial_black_hole_mass),
        np.array((initial_black_hole_angular_momentum,)),
        events=numerical_threshold,
        rtol=1.0e-10,
        atol=1.0e-12,
        max_step=0.01 * initial_black_hole_mass,
    )
    if not numerical_solution.success:
        raise RuntimeError(
            f"fixed-frequency integration failed: {numerical_solution.message}"
        )
    if len(numerical_solution.t_events[0]) != 1:
        raise RuntimeError("numerical fixed-frequency threshold was not reached")

    numerical_final_black_hole_mass = float(
        numerical_solution.t_events[0][0]
    )
    numerical_final_black_hole_angular_momentum = float(
        numerical_solution.y_events[0][0][0]
    )
    numerical_final_dimensionless_spin = (
        numerical_final_black_hole_angular_momentum
        / numerical_final_black_hole_mass**2
    )
    numerical_threshold_error = (
        fixed_frequency
        - azimuthal_number
        * horizon_frequency(
            numerical_final_black_hole_mass,
            numerical_final_dimensionless_spin
        )
    )

    return FixedFrequencyValidation(
        fixed_frequency=fixed_frequency,
        numerical_final_black_hole_mass=numerical_final_black_hole_mass,
        analytic_final_black_hole_mass=analytic_final_black_hole_mass,
        black_hole_mass_error=(
            numerical_final_black_hole_mass
            - analytic_final_black_hole_mass
        ),
        angular_momentum_error=(
            numerical_final_black_hole_angular_momentum
            - analytic_final_black_hole_angular_momentum
        ),
        numerical_threshold_error=numerical_threshold_error,
    )


def solve_backreaction_endpoint(
        alpha=0.42, initial_dimensionless_spin=0.99,
        initial_cloud_mass_fraction=1.0e-12, initial_black_hole_mass=1.0,
        ell=1, azimuthal_number=1, principal_quantum_number=2):
    """Solve the lossless single-level backreaction saturation endpoint."""
    if not 0.0 < initial_dimensionless_spin < 1.0:
        raise ValueError("initial_dimensionless_spin must lie in (0, 1)")
    if initial_cloud_mass_fraction <= 0.0:
        raise ValueError("initial_cloud_mass_fraction must be positive")
    if abs(azimuthal_number) > ell:
        raise ValueError("azimuthal_number must satisfy |m| <= ell")

    boson_mass = alpha / initial_black_hole_mass
    initial_black_hole_angular_momentum = (
        initial_dimensionless_spin * initial_black_hole_mass**2
    )
    initial_frequency = physical_mode_frequency(
        initial_black_hole_mass, initial_black_hole_angular_momentum,
        boson_mass, ell, azimuthal_number, principal_quantum_number
    )
    initial_frequency_real = initial_frequency.real
    initial_cloud_energy = (
        initial_cloud_mass_fraction * initial_black_hole_mass
    )
    initial_scaled_occupation = initial_cloud_energy / initial_frequency_real
    total_energy = initial_black_hole_mass + initial_cloud_energy
    total_angular_momentum = (
        initial_black_hole_angular_momentum
        + azimuthal_number * initial_scaled_occupation
    )

    initial_horizon_frequency = horizon_frequency(
        initial_black_hole_mass, initial_dimensionless_spin
    )
    if initial_frequency_real >= (
            azimuthal_number * initial_horizon_frequency):
        raise ValueError("the requested initial mode is not superradiant")

    def angular_momentum_derivative(black_hole_mass,
                                    angular_momentum_state):
        black_hole_angular_momentum = angular_momentum_state[0]
        frequency = physical_mode_frequency(
            black_hole_mass, black_hole_angular_momentum, boson_mass,
            ell, azimuthal_number, principal_quantum_number
        ).real
        return np.array((azimuthal_number / frequency,))

    def superradiant_threshold(black_hole_mass,
                               angular_momentum_state):
        black_hole_angular_momentum = angular_momentum_state[0]
        dimensionless_spin = (
            black_hole_angular_momentum / black_hole_mass**2
        )
        frequency = physical_mode_frequency(
            black_hole_mass, black_hole_angular_momentum, boson_mass,
            ell, azimuthal_number, principal_quantum_number
        ).real
        return (
            frequency
            - azimuthal_number
            * horizon_frequency(black_hole_mass, dimensionless_spin)
        )

    superradiant_threshold.terminal = True
    superradiant_threshold.direction = 1.0
    minimum_black_hole_mass = 0.5 * initial_black_hole_mass
    solution = solve_ivp(
        angular_momentum_derivative,
        (initial_black_hole_mass, minimum_black_hole_mass),
        np.array((initial_black_hole_angular_momentum,)),
        events=superradiant_threshold,
        rtol=1.0e-9,
        atol=1.0e-11,
        max_step=0.01 * initial_black_hole_mass,
    )
    if not solution.success:
        raise RuntimeError(f"endpoint integration failed: {solution.message}")
    if len(solution.t_events[0]) != 1:
        raise RuntimeError("superradiant threshold was not reached")

    final_black_hole_mass = float(solution.t_events[0][0])
    final_black_hole_angular_momentum = float(solution.y_events[0][0][0])
    final_dimensionless_spin = (
        final_black_hole_angular_momentum / final_black_hole_mass**2
    )
    final_frequency = physical_mode_frequency(
        final_black_hole_mass, final_black_hole_angular_momentum,
        boson_mass, ell, azimuthal_number, principal_quantum_number
    ).real
    final_horizon_frequency = horizon_frequency(
        final_black_hole_mass, final_dimensionless_spin
    )
    final_cloud_energy = (
        initial_cloud_energy
        + initial_black_hole_mass - final_black_hole_mass
    )
    final_scaled_occupation = (
        initial_scaled_occupation
        + (
            initial_black_hole_angular_momentum
            - final_black_hole_angular_momentum
        ) / azimuthal_number
    )

    final_total_energy = final_black_hole_mass + final_cloud_energy
    final_total_angular_momentum = (
        final_black_hole_angular_momentum
        + azimuthal_number * final_scaled_occupation
    )
    energy_conservation_error = final_total_energy - total_energy
    angular_momentum_conservation_error = (
        final_total_angular_momentum - total_angular_momentum
    )
    threshold_error = (
        final_frequency
        - azimuthal_number * final_horizon_frequency
    )

    initial_area = horizon_area(
        initial_black_hole_mass, initial_dimensionless_spin
    )
    final_area = horizon_area(
        final_black_hole_mass, final_dimensionless_spin
    )
    fractional_horizon_area_change = final_area / initial_area - 1.0

    maximum_residual = max(abs(energy_conservation_error),
                           abs(angular_momentum_conservation_error),
                           abs(threshold_error))
    if maximum_residual > 1.0e-9:
        raise RuntimeError(
            f"endpoint residual {maximum_residual:.3e} exceeds tolerance"
        )
    if fractional_horizon_area_change < -1.0e-10:
        raise RuntimeError("Kerr horizon area decreased during extraction")

    return BackreactionEndpoint(
        initial_black_hole_mass=initial_black_hole_mass,
        initial_dimensionless_spin=initial_dimensionless_spin,
        initial_cloud_mass_fraction=initial_cloud_mass_fraction,
        initial_scaled_occupation=initial_scaled_occupation,
        final_black_hole_mass=final_black_hole_mass,
        final_dimensionless_spin=final_dimensionless_spin,
        final_cloud_mass_fraction=(
            final_cloud_energy / final_black_hole_mass
        ),
        final_scaled_occupation=final_scaled_occupation,
        final_frequency=final_frequency,
        final_horizon_frequency=final_horizon_frequency,
        fractional_horizon_area_change=fractional_horizon_area_change,
        energy_conservation_error=energy_conservation_error,
        angular_momentum_conservation_error=(
            angular_momentum_conservation_error
        ),
        threshold_error=threshold_error,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=0.42)
    parser.add_argument("--spin", type=float, default=0.99)
    parser.add_argument("--initial-cloud-fraction", type=float, default=1.0e-12)
    parser.add_argument("--ell", type=int, default=1)
    parser.add_argument("--m", type=int, default=1)
    parser.add_argument("--principal-n", type=int, default=2)
    parser.add_argument(
        "--validate-fixed-frequency", action="store_true",
        help="compare the numerical trajectory with the exact constant-omega control"
    )
    args = parser.parse_args()

    if args.validate_fixed_frequency:
        validation = validate_fixed_frequency_trajectory(
            alpha=args.alpha,
            initial_dimensionless_spin=args.spin,
            ell=args.ell,
            azimuthal_number=args.m,
            principal_quantum_number=args.principal_n,
        )
        print("Constant-frequency exact-answer validation")
        print(f"fixed omega_R={validation.fixed_frequency:.10f}")
        print(f"endpoint M_BH: numerical="
              f"{validation.numerical_final_black_hole_mass:.12f}  "
              f"analytic={validation.analytic_final_black_hole_mass:.12f}  "
              f"error={validation.black_hole_mass_error:+.3e}")
        print(f"endpoint errors: dJ={validation.angular_momentum_error:+.3e}  "
              f"threshold={validation.numerical_threshold_error:+.3e}")
        return

    endpoint = solve_backreaction_endpoint(
        alpha=args.alpha,
        initial_dimensionless_spin=args.spin,
        initial_cloud_mass_fraction=args.initial_cloud_fraction,
        ell=args.ell,
        azimuthal_number=args.m,
        principal_quantum_number=args.principal_n,
    )

    print("Backreaction-only, lossless single-level endpoint")
    print(f"initial: M_BH={endpoint.initial_black_hole_mass:.10f}  "
          f"a*={endpoint.initial_dimensionless_spin:.10f}  "
          f"M_cloud/M_BH={endpoint.initial_cloud_mass_fraction:.4e}")
    print(f"final:   M_BH={endpoint.final_black_hole_mass:.10f}  "
          f"a*={endpoint.final_dimensionless_spin:.10f}  "
          f"M_cloud/M_BH={endpoint.final_cloud_mass_fraction:.6e}")
    print(f"threshold: omega_R={endpoint.final_frequency:.10f}  "
          f"m*Omega_H={args.m * endpoint.final_horizon_frequency:.10f}  "
          f"residual={endpoint.threshold_error:+.3e}")
    print(f"occupation: N0={endpoint.initial_scaled_occupation:.6e}  "
          f"Nf={endpoint.final_scaled_occupation:.6e}")
    print(f"horizon area change: "
          f"{endpoint.fractional_horizon_area_change:+.6e}")
    print(f"conservation: dE={endpoint.energy_conservation_error:+.3e}  "
          f"dJ={endpoint.angular_momentum_conservation_error:+.3e}")


if __name__ == "__main__":
    main()
