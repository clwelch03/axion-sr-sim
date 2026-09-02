#!/usr/bin/env python
"""Profile-level reduction gate for the multimode catalog.

The future multimode driver must reduce to the existing radial1d 211 seed when
its active mode set is ``[211]``.  This test checks the first part of that
contract: the catalog/wrapper and the established ``eigenmode_seed`` produce
the same normalized reduced radial field and normal momentum.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_RADIAL1D = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _RADIAL1D)
sys.path.insert(0, _HERE)

from radial1d import build_grid, eigenmode_seed
from mode_catalog import MODE_211, leaver_frequency, reduction_configuration
from radial_modes import build_radial_mode


def run_reduction_test(alpha=0.30, dimensionless_spin=0.90,
                       xmin=-40.0, xmax=120.0, grid_points=600):
    """Compare the multimode 211 profile with the established radial1d seed."""
    active_modes = reduction_configuration()
    if active_modes != (MODE_211,):
        raise AssertionError("single-mode reduction must activate only 211")

    x, radius = build_grid(
        dimensionless_spin, xmin, xmax, grid_points
    )
    complex_frequency = leaver_frequency(
        MODE_211, alpha, dimensionless_spin
    )
    wrapped_mode = build_radial_mode(
        MODE_211, alpha, dimensionless_spin, radius
    )
    established_reduced_profile, established_momentum_profile = eigenmode_seed(
        x, radius, alpha, dimensionless_spin, complex_frequency,
        l=MODE_211.ell, m=MODE_211.azimuthal_number, amplitude=1.0
    )

    profile_error = np.max(np.abs(
        wrapped_mode.reduced_radial_profile - established_reduced_profile
    ))
    momentum_error = np.max(np.abs(
        wrapped_mode.reduced_momentum_profile - established_momentum_profile
    ))
    if profile_error > 1.0e-12 or momentum_error > 1.0e-12:
        raise AssertionError(
            f"211 reduction mismatch: profile={profile_error:.3e}, "
            f"momentum={momentum_error:.3e}"
        )
    return profile_error, momentum_error


if __name__ == "__main__":
    profile_error, momentum_error = run_reduction_test()
    print("PASS: active mode set = [211]")
    print(f"PASS: max reduced-profile error = {profile_error:.3e}")
    print(f"PASS: max reduced-momentum error = {momentum_error:.3e}")
