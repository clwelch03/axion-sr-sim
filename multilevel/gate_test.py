#!/usr/bin/env python
"""Run all six multilevel validation gates and report a single summary.

See README.md "Required validation gates" for what each gate checks.
"""

import os
import sys
import traceback


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from mode_catalog import MODE_211, MODE_322, leaver_frequency
from angular_selection_test import run_selection_tests
from reduction_test import run_reduction_test
from evolution_reduction_test import run_evolution_reduction_test
from sign_and_hermiticity_test import run_bound_only_test, run_sign_test
from minimal_transition_test import run_minimal_transition_test


def run_label_gate():
    """Gate [1]: catalogued Leaver frequencies exist and are finite for the
    211/322 sector labels used throughout the rest of the gate suite."""
    for mode_spec in (MODE_211, MODE_322):
        frequency = leaver_frequency(mode_spec, 0.30, 0.99)
        if not (frequency == frequency):  # NaN check without numpy import
            raise AssertionError(f"{mode_spec.label} Leaver frequency is NaN")
    return True


GATES = [
    ("[1] label", run_label_gate),
    ("[2] selection", run_selection_tests),
    ("[3] bound-only", run_bound_only_test),
    ("[4] sign", run_sign_test),
    ("[5] reduction (profile)", run_reduction_test),
    ("[5] reduction (evolution)", run_evolution_reduction_test),
    ("[6] minimal transition", run_minimal_transition_test),
]


def main():
    results = []
    for name, gate in GATES:
        try:
            gate()
            print(f"PASS {name}")
            results.append(True)
        except Exception:
            print(f"FAIL {name}")
            traceback.print_exc()
            results.append(False)
    passed = sum(results)
    total = len(results)
    print(f"\n{passed}/{total} gates passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
