#!/usr/bin/env python
"""Plot the coherent |211>+|100> -> |322> population transfer demo, with the
forbidden-sector control (no |100> mediator, so |322> cannot be sourced) for
comparison.

This is the same reduced Hamiltonian toy model as
``minimal_transition_test.py``'s gate [6]: reversible, bound-only, no horizon
or continuum sink yet, and the radial overlaps use peak (not energy)
normalization.  It shows that population CAN move between hydrogenic levels
via the quartic self-interaction when the angular selection rule allows it --
not yet a calibrated transfer rate.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from transition_demo import run as run_transition


def main():
    allowed = run_transition(
        initial_amplitudes=(0.8 + 0.0j, 0.03 + 0.0j, 0.0 + 0.0j),
        final_time=300.0,
    )
    forbidden = run_transition(
        initial_amplitudes=(0.8 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j),
        final_time=300.0,
    )
    labels = [mode.label for mode in allowed["active_modes"]]
    print("active modes:", labels)
    print("max populations (allowed, 100 seeded):",
          dict(zip(labels, allowed["populations"].max(axis=0))))
    print("max populations (forbidden, no 100):",
          dict(zip(labels, forbidden["populations"].max(axis=0))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"211": "C0", "100": "C1", "322": "C2"}

    for index, label in enumerate(labels):
        ax[0].plot(allowed["t"], allowed["populations"][:, index],
                   color=colors.get(label, "k"), lw=2, label=f"|{label}>")
    ax[0].set_title(r"$211+100 \to 322$: full population")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"population $|A_j|^2$")
    ax[0].legend()

    for index, label in enumerate(labels):
        if label == "211":
            continue
        ax[1].plot(allowed["t"], allowed["populations"][:, index],
                   color=colors.get(label, "k"), lw=2, label=f"|{label}>, allowed")
    ax[1].plot(forbidden["t"], forbidden["populations"][:, 2],
               color="0.5", ls="--", lw=1.5, label="|322>, forbidden control")
    ax[1].set_yscale("log")
    ax[1].set_title("zoom: mediator + transferred population\n(forbidden control stays at machine zero)")
    ax[1].set_xlabel("t")
    ax[1].legend()

    fig.suptitle(
        "coherent bound-bound level transfer via the quartic self-interaction\n"
        "(reduced amplitude model, no horizon/continuum sink yet)"
    )
    fig.tight_layout()
    output_path = os.path.join(_HERE, "transition_demo.png")
    fig.savefig(output_path, dpi=130)
    print("saved", output_path)


if __name__ == "__main__":
    main()
