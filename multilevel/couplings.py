#!/usr/bin/env python
"""Mode-indexed angular couplings for the multilevel scaffold.

The radial mode labels live in :mod:`mode_catalog`, while the angular
projection is implemented in ``nonlinear.angular``.  This adapter keeps the
two concerns separate and returns a tensor indexed by the explicit active
``ModeSpec`` tuple.  Its convention is

    F_j = sum_{a,b,c} Q[j,a,b,c] A_a A_b^* A_c,

with ``Q[j,a,b,c]`` the angular integral ``<Y_j^* Y_a Y_b^* Y_c>``.  Entries
forbidden by azimuthal selection are exactly zero.  The radial-weighted
variant multiplies this by

    I[j,a,b,c] = integral dx u_j* u_a u_b* u_c / r(x)^2,

where ``u=r R`` is the reduced radial1d field.  The ``1/r^2`` factor comes
from projecting the cubic ``|Phi|^2 Phi`` source with ``Phi=u/r``.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_NONLINEAR = os.path.abspath(os.path.join(_HERE, "..", "nonlinear"))
_RADIAL1D = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _NONLINEAR)
sys.path.insert(0, _RADIAL1D)

from angular import quartic_overlap_general
from mode_catalog import ModeSpec, active_mode_set
from radial_modes import build_radial_mode


def quartic_coupling_tensor(active_modes):
    """Return ``(active_modes, Q)`` for an explicit active mode set.

    The tensor is deliberately angular-only.  Use
    :func:`radial_weighted_coupling_tensor` when the Leaver radial profiles
    and the reduced-equation ``1/r^2`` factor should be included.
    """
    active_modes = active_mode_set(*tuple(active_modes))
    if not all(isinstance(mode_spec, ModeSpec) for mode_spec in active_modes):
        raise TypeError("active_modes must contain only ModeSpec instances")

    mode_count = len(active_modes)
    coupling_tensor = np.zeros(
        (mode_count, mode_count, mode_count, mode_count), dtype=float
    )
    for output_index, output_mode in enumerate(active_modes):
        for first_index, first_mode in enumerate(active_modes):
            for conjugate_index, conjugate_mode in enumerate(active_modes):
                for third_index, third_mode in enumerate(active_modes):
                    coupling_tensor[output_index, first_index,
                                    conjugate_index, third_index] = (
                        quartic_overlap_general(
                            output_mode.ell,
                            output_mode.azimuthal_number,
                            first_mode.ell,
                            first_mode.azimuthal_number,
                            conjugate_mode.ell,
                            conjugate_mode.azimuthal_number,
                            third_mode.ell,
                            third_mode.azimuthal_number,
                        )
                    )
    return active_modes, coupling_tensor


def cubic_force(amplitudes, coupling_tensor):
    """Evaluate ``F_j = Q[j,a,b,c] A_a A_b* A_c``."""
    amplitudes = np.asarray(amplitudes, dtype=complex)
    if coupling_tensor.shape != (len(amplitudes),) * 4:
        raise ValueError("coupling tensor shape does not match amplitudes")
    return np.einsum(
        "jabc,a,b,c->j", coupling_tensor, amplitudes,
        np.conjugate(amplitudes), amplitudes,
    )


def radial_weighted_coupling_tensor(active_modes, alpha, dimensionless_spin,
                                    x, radius, series_terms=400):
    """Return ``(modes, radial_modes, Q_full)`` including radial overlaps.

    The returned profiles use the existing peak normalization
    ``max(abs(u))=1``.  That keeps this tensor comparable to the current
    radial1d amplitude convention, but it is not yet an energy-normalized
    production coupling.
    """
    active_modes, angular_tensor = quartic_coupling_tensor(active_modes)
    x_array = np.asarray(x, dtype=float)
    radius_array = np.asarray(radius, dtype=float)
    if (x_array.ndim != 1 or radius_array.ndim != 1
            or x_array.shape != radius_array.shape
            or x_array.size < 2):
        raise ValueError("x and radius must be matching one-dimensional grids")
    if (not np.all(np.isfinite(x_array))
            or not np.all(np.isfinite(radius_array))
            or np.any(np.diff(x_array) <= 0.0)
            or np.any(radius_array <= 0.0)):
        raise ValueError("x and radius must be finite, ordered, and exterior")

    radial_modes = tuple(
        build_radial_mode(
            mode_spec, alpha, dimensionless_spin, radius_array,
            series_terms=series_terms,
        )
        for mode_spec in active_modes
    )
    reduced_profiles = np.asarray([
        radial_mode.reduced_radial_profile for radial_mode in radial_modes
    ])
    radial_tensor = np.zeros_like(angular_tensor, dtype=complex)
    radius_factor = radius_array ** 2
    for output_index in range(len(active_modes)):
        for first_index in range(len(active_modes)):
            for conjugate_index in range(len(active_modes)):
                for third_index in range(len(active_modes)):
                    if angular_tensor[output_index, first_index,
                                      conjugate_index, third_index] == 0.0:
                        continue
                    integrand = (
                        np.conjugate(reduced_profiles[output_index])
                        * reduced_profiles[first_index]
                        * np.conjugate(reduced_profiles[conjugate_index])
                        * reduced_profiles[third_index]
                        / radius_factor
                    )
                    radial_tensor[output_index, first_index,
                                  conjugate_index, third_index] = np.trapezoid(
                        integrand, x_array
                    )
    return active_modes, radial_modes, angular_tensor * radial_tensor


def _selftest():
    from mode_catalog import MODE_211, MODE_322, ModeSpec

    zero_mode = ModeSpec(1, 0, 0)
    modes, tensor = quartic_coupling_tensor((MODE_211, zero_mode, MODE_322))
    labels = [mode.label for mode in modes]
    allowed = tensor[2, 0, 1, 0]
    forbidden = tensor[2, 0, 0, 0]
    print("active modes:", labels)
    print("Q[322,211,100,211] =", allowed, "(want nonzero)")
    print("Q[322,211,211,211] =", forbidden, "(want 0)")
    if abs(allowed) == 0.0 or forbidden != 0.0:
        raise AssertionError("multilevel selection-rule adapter failed")

    from radial1d import build_grid
    x, radius = build_grid(0.99, -20.0, 60.0, 160)
    _, radial_modes, full_tensor = radial_weighted_coupling_tensor(
        modes, 0.30, 0.99, x, radius, series_terms=120
    )
    radial_allowed = full_tensor[2, 0, 1, 0]
    print("radial-weighted Q[322,211,100,211] =", radial_allowed)
    if not np.isfinite(radial_allowed):
        raise AssertionError("radial overlap is non-finite")
    if not all(np.max(np.abs(mode.reduced_radial_profile)) <= 1.0 + 1e-12
               for mode in radial_modes):
        raise AssertionError("radial profile normalization failed")


if __name__ == "__main__":
    _selftest()
