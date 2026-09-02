#!/usr/bin/env python
"""Full-evolution 211 reduction gate for the multimode linear driver."""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_RADIAL1D = os.path.abspath(os.path.join(_HERE, ".."))
_NONLINEAR = os.path.join(_RADIAL1D, "nonlinear")
sys.path.insert(0, _HERE)
sys.path.insert(0, _NONLINEAR)

from evolve_tower import run as run_existing_tower
from linear_driver import run as run_multilevel_linear
from mode_catalog import MODE_211, reduction_configuration


def run_evolution_reduction_test():
    """Compare the independent-channel driver with the established l=1 tower."""
    parameters = dict(
        alpha=0.30,
        dimensionless_spin=0.90,
        xmin=-35.0,
        xmax=90.0,
        grid_points=400,
        final_time=60.0,
        eps_match=0.02,
        pml_peak_damping=1.0,
        recorded_points=120,
    )
    multimode_output = run_multilevel_linear(
        reduction_configuration(), **parameters
    )
    tower_output = run_existing_tower(
        alpha=parameters["alpha"],
        a=parameters["dimensionless_spin"],
        l_max=1,
        xmin=parameters["xmin"],
        xmax=parameters["xmax"],
        N=parameters["grid_points"],
        T=parameters["final_time"],
        eps_match=parameters["eps_match"],
        gmax=parameters["pml_peak_damping"],
        nrec=parameters["recorded_points"],
        lambda_quartic=0.0,
        init_amplitude=1.0,
        seed="eigenmode",
    )

    if not np.array_equal(multimode_output["t"], tower_output["t"]):
        raise AssertionError("multimode and tower recording times differ")
    field_error = np.max(np.abs(
        multimode_output["probe_field"] - tower_output["psi"]
    ))
    norm_error = np.max(np.abs(
        multimode_output["cloud_norm"] - tower_output["cloud_norm"]
    ))
    if field_error > 1.0e-11 or norm_error > 1.0e-11:
        raise AssertionError(
            f"211 evolution mismatch: field={field_error:.3e}, "
            f"cloud_norm={norm_error:.3e}"
        )
    return field_error, norm_error


if __name__ == "__main__":
    field_error, norm_error = run_evolution_reduction_test()
    print("PASS: multimode active mode set = [211]")
    print(f"PASS: max probe-field error = {field_error:.3e}")
    print(f"PASS: max cloud-norm error = {norm_error:.3e}")
