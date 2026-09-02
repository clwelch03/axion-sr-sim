#!/usr/bin/env python
"""Plot the effect of a non-Hermitian horizon/continuum sink on the |322>
channel of the reduced |211>+|100> -> |322> transition demo.

Same reduced amplitude model and parameters as ``plot_transition_demo.py``'s
allowed run, but here we add an explicit sink ``loss_rates=[0, 0, Gamma_322]``
on the |322> channel only (a stand-in for the horizon/continuum flux that the
Hermitian bound-only model omits). Comparing against the reversible
(loss_rates=None) reference shows the qualitative change from a Rabi-like
back-and-forth sloshing of population into a secular one-way drain of
|211>'s occupation: once population is pumped into |322>, the sink removes it
before it can coherently transfer back, so each Rabi pulse leaves |211> a
little lower than before.

``Gamma_322`` is picked relative to the Rabi angular frequency set by the
~40M pulse period visible in transition_demo.png, omega_Rabi ~ 2*pi/40 ~
0.157. A first scan at Gamma_322 ~ omega_Rabi (and a few times larger) only
drained |211> by ~1-2% over t=300 -- the |322> occupation is tiny (peak
~0.003) so a weak sink barely dents |211>. A follow-up scan up to
Gamma_322 ~ 10*omega_Rabi showed the drain fraction actually peaks around
Gamma_322 ~ 1.0-1.5 (~6-10x omega_Rabi) at ~2.6-2.8%, then falls off for
larger rates: once the sink is much faster than the coupling can refill
|322>, it starts adiabatically eliminating (Zeno-suppressing) the channel
before population can build up in it at all, so very strong sinks transfer
*less* net population out of |211>, not more. The candidates below bracket
that peak; we keep whichever shows the clearest monotonic (secular) decline
of |211> against the reversible reference's ~40M oscillation.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from transition_demo import run as run_transition

COMMON_KWARGS = dict(
    initial_amplitudes=(0.8 + 0.0j, 0.03 + 0.0j, 0.0 + 0.0j),
    final_time=300.0,
)

# omega_Rabi ~ 2*pi / (~40M pulse period) ~ 0.157 (see transition_demo.png).
# The |211> drain fraction was scanned separately over Gamma_322 in
# [0.157, 10] and peaks around Gamma_322 ~ 1.0-1.5 (~6-10x omega_Rabi); these
# candidates bracket that peak.
OMEGA_RABI = 2.0 * np.pi / 40.0
SINK_CANDIDATES = {
    "at peak (~6.4x omega_Rabi)": 1.0,
    "above peak (~9.5x omega_Rabi)": 1.5,
    "well above peak (~19x omega_Rabi)": 3.0,
}


def drain_metric(reference, sunk):
    """Fractional drop of |211> population at t=300 relative to the
    reversible reference, plus the number of resolved |322> pulses (local
    maxima) so we can check the sink isn't erasing pulse structure."""
    population_211_reference = reference["populations"][-1, 0]
    population_211_sunk = sunk["populations"][-1, 0]
    fractional_drain = 1.0 - population_211_sunk / population_211_reference
    population_322_sunk = sunk["populations"][:, 2]
    is_local_max = (
        (population_322_sunk[1:-1] > population_322_sunk[:-2])
        & (population_322_sunk[1:-1] > population_322_sunk[2:])
        & (population_322_sunk[1:-1] > 1e-6)
    )
    pulse_count = int(np.sum(is_local_max))
    return fractional_drain, pulse_count


def main():
    print("running reversible reference (loss_rates=None) ...")
    reference = run_transition(loss_rates=None, **COMMON_KWARGS)
    labels = [mode.label for mode in reference["active_modes"]]
    print("active modes:", labels)

    candidate_results = {}
    for name, gamma_322 in SINK_CANDIDATES.items():
        print(f"running sink candidate {name!r}: Gamma_322={gamma_322:.4f} ...")
        result = run_transition(
            loss_rates=[0.0, 0.0, gamma_322], **COMMON_KWARGS
        )
        fractional_drain, pulse_count = drain_metric(reference, result)
        candidate_results[name] = (gamma_322, result, fractional_drain, pulse_count)
        print(f"    Gamma_322={gamma_322:.4f}: |211> fractional drain at t=300 "
              f"= {fractional_drain:.3f}, resolved |322> pulses = {pulse_count}")

    # Pick the candidate with the largest visible drain that still leaves at
    # least a couple of resolved |322> pulses (so the panel shows damped
    # pulses, not just an instant wipeout on step one).
    eligible = {
        name: vals for name, vals in candidate_results.items() if vals[3] >= 2
    }
    if not eligible:
        eligible = candidate_results
    chosen_name = max(eligible, key=lambda name: eligible[name][2])
    gamma_322, sunk, fractional_drain, pulse_count = candidate_results[chosen_name]
    print(f"chosen sink candidate: {chosen_name!r}, Gamma_322={gamma_322:.4f}, "
          f"fractional |211> drain={fractional_drain:.3f}, pulses={pulse_count}")

    print("max populations (reference, no sink):",
          dict(zip(labels, reference["populations"].max(axis=0))))
    print("max populations (with sink):",
          dict(zip(labels, sunk["populations"].max(axis=0))))
    print("final populations (reference, no sink):",
          dict(zip(labels, reference["populations"][-1])))
    print("final populations (with sink):",
          dict(zip(labels, sunk["populations"][-1])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"211": "C0", "100": "C1", "322": "C2"}
    index_211 = labels.index("211")
    index_322 = labels.index("322")

    ax[0].plot(reference["t"], reference["populations"][:, index_211],
               color=colors["211"], lw=2, ls="--", label="|211>, no sink")
    ax[0].plot(sunk["t"], sunk["populations"][:, index_211],
               color=colors["211"], lw=2, label="|211>, with |322> sink")
    ax[0].set_title(r"$|211\rangle$ population: reversible vs. sunk $|322\rangle$")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"population $|A_{211}|^2$")
    ax[0].legend()

    ax[1].plot(reference["t"], reference["populations"][:, index_322],
               color=colors["322"], lw=1.5, ls="--", label="|322>, no sink")
    ax[1].plot(sunk["t"], sunk["populations"][:, index_322],
               color=colors["322"], lw=2, label="|322>, with sink")
    ax[1].set_yscale("log")
    ax[1].set_title(r"$|322\rangle$ population (log): periodic pulses vs. damped drain")
    ax[1].set_xlabel("t")
    ax[1].legend()

    fig.suptitle(
        "non-Hermitian horizon/continuum sink on |322> turns Rabi-like\n"
        "sloshing into secular drain of |211> "
        rf"($\Gamma_{{322}}$={gamma_322:.3f}, {chosen_name})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    output_path = os.path.join(_HERE, "transition_sink_demo.png")
    fig.savefig(output_path, dpi=130)
    print("saved", output_path)


if __name__ == "__main__":
    main()
