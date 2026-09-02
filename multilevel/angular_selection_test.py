#!/usr/bin/env python
"""Selection-rule tests for arbitrary signed-m quartic overlaps."""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "nonlinear"))

import angular


def run_selection_tests():
    common_m_overlap = angular.quartic_overlap(1, 1, 1, 1, m=1)
    general_m_overlap = angular.quartic_overlap_general(
        1, 1, 1, 1, 1, 1, 1, 1
    )
    np.testing.assert_allclose(common_m_overlap, general_m_overlap, rtol=0.0, atol=1e-15)

    # 211 + 211 with an m=0 conjugate factor can source an ell=2,m=2 output.
    allowed_overlap = angular.quartic_overlap_general(
        2, 2, 1, 1, 0, 0, 1, 1
    )
    if abs(allowed_overlap) <= 1.0e-12:
        raise AssertionError("allowed 211/100/211 -> 322 overlap vanished")

    forbidden_overlap = angular.quartic_overlap_general(
        2, 1, 1, 1, 0, 0, 1, 1
    )
    if forbidden_overlap != 0.0:
        raise AssertionError("forbidden m-selection channel was not zero")

    # The signed-m theta factors preserve spherical-harmonic normalization.
    positive_norm = 2.0 * np.pi * np.sum(
        angular._W * angular._Lambda(1, 1, angular._U) ** 2
    )
    negative_norm = 2.0 * np.pi * np.sum(
        angular._W * angular._Lambda(1, -1, angular._U) ** 2
    )
    np.testing.assert_allclose(positive_norm, 1.0, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(negative_norm, 1.0, rtol=0.0, atol=1e-14)
    return common_m_overlap, allowed_overlap


if __name__ == "__main__":
    common_m_overlap, allowed_overlap = run_selection_tests()
    print("PASS: common-m API matches generalized overlap")
    print(f"PASS: allowed 211/100/211 -> 322 overlap = {allowed_overlap:+.6e}")
    print("PASS: forbidden m-selection channel = 0")
    print("PASS: signed-m spherical-harmonic normalization")
