#!/usr/bin/env python
"""Leaver radial-mode wrapper for the multilevel scaffold.

The wrapper is intentionally a profile builder, not an evolver.  It keeps the
two radial conventions visible:

* ``physical_radial_profile`` is the separated radial function R(r) in
  ``Phi = R(r) Y_lm exp(-i omega t)``;
* ``reduced_radial_profile`` is r R(r), the radial1d field ``psi`` in
  ``Phi = psi(t,x) Y_lm / r``.

The default normalization is ``max(abs(r R)) = 1`` over the supplied grid,
matching ``radial1d.eigenmode_seed``.  A production mode-coupling calculation
must replace this with an explicitly documented energy normalization.
"""

from dataclasses import dataclass
import os
import sys

import mpmath as mp
import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_LINEAR_RATES = os.path.abspath(
    os.path.join(_HERE, "..", "..", "linear_rates")
)
sys.path.insert(0, _LINEAR_RATES)

import leaver_eigenfunction as leaver_radial

from mode_catalog import ModeSpec, leaver_frequency


@dataclass(frozen=True)
class LeaverRadialMode:
    """A normalized Leaver mode sampled on a supplied BL-radius grid."""

    mode_spec: ModeSpec
    dimensionless_mass_coupling: float
    dimensionless_spin: float
    complex_frequency: complex
    radius: np.ndarray
    physical_radial_profile: np.ndarray
    reduced_radial_profile: np.ndarray
    reduced_momentum_profile: np.ndarray


def build_radial_mode(mode_spec, dimensionless_mass_coupling,
                      dimensionless_spin, radius, series_terms=400):
    """Build one normalized Leaver radial mode on a BL-radius grid.

    Parameters
    ----------
    mode_spec : ModeSpec
        Full `(principal_n, ell, m)` channel identity.
    dimensionless_mass_coupling : float
        `M * mu`, with the mode and radius convention using `M=1`.
    dimensionless_spin : float
        Kerr `chi = a/M` in `(0, 1)`.
    radius : array_like
        Strictly exterior Boyer-Lindquist radii, larger than `r_plus`.
    series_terms : int
        Number of Leaver radial-series terms retained before the tail.
    """
    radius_array = np.asarray(radius, dtype=float)
    if radius_array.ndim != 1 or radius_array.size < 2:
        raise ValueError("radius must be a one-dimensional grid with at least two points")
    if np.any(~np.isfinite(radius_array)) or np.any(np.diff(radius_array) <= 0.0):
        raise ValueError("radius must be finite and strictly increasing")
    if not 0.0 < dimensionless_spin < 1.0:
        raise ValueError("dimensionless_spin must lie in (0, 1)")
    horizon_radius = 1.0 + np.sqrt(1.0 - dimensionless_spin**2)
    if radius_array[0] <= horizon_radius:
        raise ValueError("radius grid must lie strictly outside the outer horizon")
    if series_terms < 10:
        raise ValueError("series_terms is too small for a radial profile")

    complex_frequency = leaver_frequency(
        mode_spec, dimensionless_mass_coupling, dimensionless_spin
    )
    series_coefficients = leaver_radial.series_coeffs(
        mp.mpc(complex_frequency), dimensionless_mass_coupling,
        dimensionless_spin, mode_spec.ell, mode_spec.azimuthal_number,
        series_terms
    )
    physical_radial_profile = np.array([
        complex(leaver_radial.R_of_r(
            mp.mpf(radius_value), mp.mpc(complex_frequency),
            dimensionless_mass_coupling, dimensionless_spin,
            series_coefficients, m=mode_spec.azimuthal_number
        ))
        for radius_value in radius_array
    ], dtype=complex)

    reduced_radial_profile = radius_array * physical_radial_profile
    normalization = np.max(np.abs(reduced_radial_profile))
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("Leaver radial profile has no finite nonzero normalization")
    reduced_radial_profile = reduced_radial_profile / normalization
    physical_radial_profile = physical_radial_profile / normalization
    reduced_momentum_profile = -1j * complex_frequency * reduced_radial_profile

    return LeaverRadialMode(
        mode_spec=mode_spec,
        dimensionless_mass_coupling=float(dimensionless_mass_coupling),
        dimensionless_spin=float(dimensionless_spin),
        complex_frequency=complex_frequency,
        radius=radius_array,
        physical_radial_profile=physical_radial_profile,
        reduced_radial_profile=reduced_radial_profile,
        reduced_momentum_profile=reduced_momentum_profile,
    )


if __name__ == "__main__":
    from mode_catalog import MODE_211, MODE_322

    sample_radius = np.linspace(2.0, 40.0, 128)
    for mode_spec in (MODE_211, MODE_322):
        mode = build_radial_mode(mode_spec, 0.30, 0.99, sample_radius)
        print(
            f"{mode_spec.label}: omega={mode.complex_frequency.real:.12f}"
            f"{mode.complex_frequency.imag:+.3e}i, "
            f"max|rR|={np.max(np.abs(mode.reduced_radial_profile)):.6f}"
        )
