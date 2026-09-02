#!/usr/bin/env python
"""Minimal transition gate [6]: 211+100 -> 322 populates; a forbidden-sector
control (211 alone, no 100 mediator) leaves 322 unpopulated.

This exercises the full coherent evolution in :mod:`transition_demo`, not
just the static angular tensor (already checked by
``angular_selection_test.py`` and ``couplings.py``'s self-test).  The
distinction matters because a nonzero |100> seed is what actually opens a
term with nonzero azimuthal-selection weight in the equation of motion for
the |322> channel; with |100> absent, that channel's forcing vanishes even
though 322 sits in the active mode set with its own linear frequency.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from transition_demo import run as run_transition


def run_minimal_transition_test():
    allowed_run = run_transition(
        initial_amplitudes=(0.8 + 0.0j, 0.03 + 0.0j, 0.0 + 0.0j),
        final_time=300.0,
    )
    forbidden_run = run_transition(
        initial_amplitudes=(0.8 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j),
        final_time=300.0,
    )

    allowed_322_population = allowed_run["populations"][:, 2].max()
    forbidden_322_population = forbidden_run["populations"][:, 2].max()

    if not allowed_322_population > 1.0e-4:
        raise AssertionError(
            "211+100 -> 322 transfer failed to populate the 322 channel: "
            f"max population = {allowed_322_population:.3e}"
        )
    if not forbidden_322_population < 1.0e-12:
        raise AssertionError(
            "forbidden-sector control (211 alone) spuriously populated 322: "
            f"max population = {forbidden_322_population:.3e}"
        )
    return allowed_322_population, forbidden_322_population


if __name__ == "__main__":
    allowed_322_population, forbidden_322_population = run_minimal_transition_test()
    print(
        f"PASS: allowed 211+100 -> 322 max population = "
        f"{allowed_322_population:.3e}"
    )
    print(
        f"PASS: forbidden-sector (211 alone) max 322 population = "
        f"{forbidden_322_population:.3e}"
    )
