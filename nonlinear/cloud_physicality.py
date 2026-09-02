#!/usr/bin/env python
"""Translate the radial1d seed knob into cloud size and mass diagnostics.

The evolution uses the complex envelope

    Phi_complex = (psi / r) Y_lm exp(i m phi),

with a physical real axion represented as Re(Phi_complex).  The seed is
normalised by ``max|psi| = A``.  Consequently, A is neither the peak axion
angle nor a cloud mass.

This script reconstructs the clean Leaver seed and evaluates the exterior Kerr
Killing energy from Dolan 2012, arXiv:1212.1477, Eqs. (27)--(30).  It reports:

1. the physical-field peak and energy obtained by interpreting A directly as
   the radial amplitude of the dimensionless axion angle;
2. the axion amplitude with the same |lambda| A^2 as a toy-code run.

The second comparison matches only the magnitude of the nonlinear force.  A
cosine axion is attractive, whereas the existing Gamma(A,l_max) sweep used the
repulsive sign lambda_quartic=+0.02.
"""

import argparse
import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_RADIAL1D = os.path.dirname(_HERE)
_LINEAR_RATES = os.path.abspath(os.path.join(_RADIAL1D, "..", "linear_rates"))
sys.path.insert(0, _RADIAL1D)
sys.path.insert(0, _LINEAR_RATES)

import angular
import leaver
from radial1d import (build_grid, d_dx, eigenmode_seed, horizons,
                      tortoise_x)


def dolan_complex_energy_density(x, radius, psi, momentum, alpha, spin,
                                 ell=1, azimuthal_number=1):
    """Return Dolan's complex-field exterior Killing-energy density dE/dx.

    This is the diagonal, single-l form of arXiv:1212.1477 Eq. (29), consistent
    with the linear operator used by radial1d.  ``psi`` is the solver's reduced
    radial field; the physical radial coefficient is ``psi / r``.

    The nonlinear quartic potential energy is deliberately excluded.  This
    diagnostic establishes the seed/cloud normalization against the linear
    quasi-bound state; the production sweep used a toy repulsive RWA coupling,
    not the physical cosine-axion potential.
    """
    grid_spacing = x[1] - x[0]
    delta = radius**2 - 2.0 * radius + spin**2
    radial_factor = radius**2 + spin**2
    sigma_zero_squared = radial_factor**2 - spin**2 * delta
    cos2_diagonal = angular.cos2_coupling(ell, ell, azimuthal_number)
    time_coefficient = (
        sigma_zero_squared + delta * spin**2 * cos2_diagonal
    )

    physical_radial_field = psi / radius
    physical_radial_momentum = momentum / radius
    radial_derivative = d_dx(physical_radial_field, grid_spacing)

    twice_radial_factor_times_density = (
        delta
        * (ell * (ell + 1)
           + alpha**2 * (radius**2 + spin**2 * cos2_diagonal))
        * np.abs(physical_radial_field)**2
        + time_coefficient * np.abs(physical_radial_momentum)**2
        + radial_factor**2 * np.abs(radial_derivative)**2
        + 2.0 * spin * azimuthal_number * radial_factor
        * np.imag(np.conj(physical_radial_field) * radial_derivative)
    )
    return twice_radial_factor_times_density / (2.0 * radial_factor)


def cumulative_energy_radii(radius, energy_density, x, fractions):
    """Return radii enclosing selected fractions of positive exterior energy."""
    shell_energy = 0.5 * (energy_density[1:] + energy_density[:-1]) * np.diff(x)
    cumulative_energy = np.concatenate(([0.0], np.cumsum(shell_energy)))
    if cumulative_energy[-1] <= 0.0:
        raise ValueError("Exterior Killing energy is not positive")
    cumulative_fraction = cumulative_energy / cumulative_energy[-1]
    if np.any(np.diff(cumulative_fraction) < -1e-12):
        raise ValueError("Energy density is not monotone-positive on this domain")
    return {
        fraction: float(np.interp(fraction, cumulative_fraction, radius))
        for fraction in fractions
    }


def seed_diagnostics(alpha=0.42, spin=0.99, amplitudes=(1.0, 2.0, 4.0),
                     lambda_quartic=0.02, axion_decay_constant_gev=1.0e16,
                     planck_mass_gev=1.22089e19, xmin=-750.0, xmax=250.0,
                     grid_points=6000, eps_match=0.002,
                     outer_sponge_width=40.0):
    """Compute physicality diagnostics for the clean l=m=1 eigenmode seed."""
    ell = azimuthal_number = 1
    omega = complex(leaver.solve(alpha, spin))
    x, radius = build_grid(spin, xmin, xmax, grid_points)
    unit_psi, unit_momentum = eigenmode_seed(
        x, radius, alpha, spin, omega, ell, azimuthal_number, amplitude=1.0
    )

    horizon_radius, _ = horizons(spin)
    x_match = tortoise_x(horizon_radius + eps_match, spin)
    outer_sponge_start = xmax - outer_sponge_width
    physical_region = (x >= x_match) & (x <= outer_sponge_start)

    energy_density = dolan_complex_energy_density(
        x, radius, unit_psi, unit_momentum, alpha, spin,
        ell, azimuthal_number
    )
    complex_envelope_energy = np.trapezoid(
        energy_density[physical_region], x[physical_region]
    )
    # For phi_real = Re(phi_complex), time/azimuth averaging gives half the
    # energy of Dolan's canonically normalised complex field.
    real_field_energy_coefficient = 0.5 * complex_envelope_energy

    y11_peak = np.sqrt(3.0 / (8.0 * np.pi))
    field_envelope = np.abs(unit_psi / radius) * y11_peak
    physical_field_peak_per_amplitude = float(np.max(field_envelope[physical_region]))
    peak_index = np.flatnonzero(physical_region)[
        np.argmax(field_envelope[physical_region])
    ]

    enclosed_radii = cumulative_energy_radii(
        radius[physical_region], energy_density[physical_region],
        x[physical_region], fractions=(0.10, 0.50, 0.90, 0.99)
    )

    quartic_overlap_1111 = angular.quartic_overlap(1, 1, 1, 1, 1)
    physical_axion_lambda = -alpha**2 * quartic_overlap_1111 / 8.0
    nonlinear_strength_ratio = abs(lambda_quartic / physical_axion_lambda)
    decay_constant_ratio_squared = (
        axion_decay_constant_gev / planck_mass_gev
    )**2

    rows = []
    for code_amplitude in amplitudes:
        matched_axion_amplitude = (
            code_amplitude * np.sqrt(nonlinear_strength_ratio)
        )
        direct_mass_fraction = (
            real_field_energy_coefficient * code_amplitude**2
            * decay_constant_ratio_squared
        )
        matched_mass_fraction = (
            real_field_energy_coefficient * matched_axion_amplitude**2
            * decay_constant_ratio_squared
        )
        rows.append({
            "code_amplitude": code_amplitude,
            "code_nonlinear_strength": lambda_quartic * code_amplitude**2,
            "direct_peak_axion_angle": (
                physical_field_peak_per_amplitude * code_amplitude
            ),
            "direct_mass_fraction": direct_mass_fraction,
            "matched_axion_amplitude": matched_axion_amplitude,
            "matched_peak_axion_angle": (
                physical_field_peak_per_amplitude * matched_axion_amplitude
            ),
            "matched_mass_fraction": matched_mass_fraction,
        })

    return {
        "omega": omega,
        "physical_field_peak_per_amplitude": physical_field_peak_per_amplitude,
        "physical_field_peak_radius": float(radius[peak_index]),
        "complex_envelope_energy": complex_envelope_energy,
        "real_field_energy_coefficient": real_field_energy_coefficient,
        "enclosed_radii": enclosed_radii,
        "physical_axion_lambda": physical_axion_lambda,
        "nonlinear_strength_ratio": nonlinear_strength_ratio,
        "decay_constant_ratio_squared": decay_constant_ratio_squared,
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=0.42)
    parser.add_argument("--spin", type=float, default=0.99)
    parser.add_argument("--amplitudes", type=float, nargs="+", default=(1.0, 2.0, 4.0))
    parser.add_argument("--lambda", dest="lambda_quartic", type=float, default=0.02)
    parser.add_argument("--fa-gev", type=float, default=1.0e16)
    parser.add_argument("--planck-mass-gev", type=float, default=1.22089e19)
    args = parser.parse_args()

    result = seed_diagnostics(
        alpha=args.alpha, spin=args.spin, amplitudes=args.amplitudes,
        lambda_quartic=args.lambda_quartic,
        axion_decay_constant_gev=args.fa_gev,
        planck_mass_gev=args.planck_mass_gev,
    )

    print(f"Leaver omega = {result['omega'].real:.9f} "
          f"{result['omega'].imag:+.4e} i")
    print("Seed convention: max|psi|=A, physical axion angle=Re[(psi/r)Y_11 e^{i phi}]")
    print(f"peak |Phi|/f_a per A = {result['physical_field_peak_per_amplitude']:.8g} "
          f"at r={result['physical_field_peak_radius']:.3f} M")
    print(f"Dolan complex exterior energy / A^2 = {result['complex_envelope_energy']:.8g}")
    print(f"real-envelope averaged energy / A^2 = "
          f"{result['real_field_energy_coefficient']:.8g}")
    print("energy enclosure radii: " + ", ".join(
        f"r{int(100*fraction)}={radius:.2f} M"
        for fraction, radius in result["enclosed_radii"].items()
    ))
    print(f"cosine-axion lambda_code = {result['physical_axion_lambda']:+.8g}")
    print(f"|lambda_toy/lambda_axion| = {result['nonlinear_strength_ratio']:.6g}")
    print(f"(f_a/M_Pl)^2 = {result['decay_constant_ratio_squared']:.8g}")
    if args.lambda_quartic * result["physical_axion_lambda"] < 0.0:
        print("WARNING: toy run and cosine axion have opposite interaction signs")

    heading = (
        "A_code  lambda*A^2  direct peak  direct Mcloud/M  "
        "A_axion(|g| match)  matched peak  matched Mcloud/M"
    )
    print("\n" + heading)
    for row in result["rows"]:
        print(
            f"{row['code_amplitude']:6.2f}  "
            f"{row['code_nonlinear_strength']:10.4g}  "
            f"{row['direct_peak_axion_angle']:11.5g}  "
            f"{row['direct_mass_fraction']:15.5g}  "
            f"{row['matched_axion_amplitude']:18.5g}  "
            f"{row['matched_peak_axion_angle']:12.5g}  "
            f"{row['matched_mass_fraction']:16.5g}"
        )


if __name__ == "__main__":
    main()
