#!/usr/bin/env python
"""Bound-only Hermiticity gate [3] and horizon/continuum sign gate [4].

Both gates use the reduced three-mode (|211>, |100>, |322>) transition model
in :mod:`transition_demo`.  They are independent of whether the bound-bound
coupling tensor is Hermitian in isolation (already covered by the exact
instantaneous identity in ``nonlinear/angular.py``'s tower gate); here the
question is whether the *retained-mode* system as a whole conserves norm when
no sink is present, and whether an explicit sink of a given sign moves
occupation the right way.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from transition_demo import run as run_transition


def run_bound_only_test():
    """Gate [3]: with every loss_rate zero, total norm is conserved."""
    output = run_transition(loss_rates=None, final_time=300.0)
    norm_drift = np.max(np.abs(output["norm"] - output["norm"][0]))
    if norm_drift > 1.0e-9:
        raise AssertionError(
            f"bound-only system leaked norm: max drift = {norm_drift:.3e}"
        )
    return norm_drift


def run_sign_test():
    """Gate [4]: a positive loss rate depletes its channel; a negative one
    (gain) does not.  The cubic coupling is disabled so the two runs isolate
    the sign of the explicit sink term rather than any coherent transfer."""
    channel_count = 3  # (211, 100, 322), matching transition_demo's active set
    loss_rate = 5.0e-3
    positive_loss = np.zeros(channel_count)
    positive_loss[0] = loss_rate
    negative_loss = np.zeros(channel_count)
    negative_loss[0] = -loss_rate

    depleting_run = run_transition(
        interaction_strength=0.0, loss_rates=positive_loss, final_time=200.0
    )
    growing_run = run_transition(
        interaction_strength=0.0, loss_rates=negative_loss, final_time=200.0
    )

    initial_population = depleting_run["populations"][0, 0]
    depleted_population = depleting_run["populations"][-1, 0]
    grown_population = growing_run["populations"][-1, 0]

    if not depleted_population < initial_population:
        raise AssertionError(
            "positive loss_rate failed to deplete the sourced channel: "
            f"{initial_population:.6e} -> {depleted_population:.6e}"
        )
    if not grown_population > initial_population:
        raise AssertionError(
            "negative loss_rate (gain) failed to grow the sourced channel: "
            f"{initial_population:.6e} -> {grown_population:.6e}"
        )
    expected_decay = initial_population * np.exp(-2.0 * loss_rate * 200.0)
    decay_error = abs(depleted_population - expected_decay) / expected_decay
    if decay_error > 1.0e-6:
        raise AssertionError(
            f"loss-only channel did not follow exp(-2*loss_rate*t): "
            f"relative error = {decay_error:.3e}"
        )
    return initial_population, depleted_population, grown_population


if __name__ == "__main__":
    norm_drift = run_bound_only_test()
    print(f"PASS: bound-only norm drift = {norm_drift:.3e}")
    initial_population, depleted_population, grown_population = run_sign_test()
    print(
        "PASS: sign gate -- 211 population "
        f"{initial_population:.6e} -> depleted {depleted_population:.6e} "
        f"(+loss_rate), -> grown {grown_population:.6e} (-loss_rate)"
    )
