"""Plot the four reference curves exported by reference_profiles.jl.

Usage: PYTHONPATH=src python plot_reference_profiles.py RESULT_PREFIX FIGURE_PREFIX
The CSV contains physical radial coefficients; no normalization is applied here.
"""

import argparse
import csv
from pathlib import Path
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import plotstyle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_prefix", type=Path)
    parser.add_argument("figure_prefix", type=Path)
    args = parser.parse_args()

    with args.result_prefix.with_suffix(".toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    normalization = metadata["normalization"]
    if normalization not in ("fixed_outer_coefficients", "equal_cloud_energies"):
        raise ValueError("Unrecognized energy/amplitude normalization.")

    # Preserve both matching-point samples from the CSV. A jump is not smoothed.
    curves = {}
    with args.result_prefix.with_suffix(".csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            curves.setdefault(row["curve"], []).append(
                (float(row["r_star_over_M"]), float(row["R_abs"]))
            )

    styles = (
        ("omiya_primary", r"Omiya $l=m=1$", plotstyle.PALETTE[4], "-"),
        ("linear_primary", r"Linear $l=m=1$", plotstyle.PALETTE[3], ":"),
        ("omiya_secondary", r"Omiya $l=m=2$", plotstyle.PALETTE[1], "-"),
        ("linear_secondary", r"Linear $l=m=2$", plotstyle.PALETTE[0], "--"),
    )
    if set(curves) != {name for name, *_ in styles}:
        raise ValueError("Expected exactly the four reference curves.")

    fig, ax = plt.subplots(figsize=(9, 6))
    for name, label, color, linestyle in styles:
        samples = np.asarray(curves[name])
        if not np.isfinite(samples).all() or np.any(np.diff(samples[:, 0]) < 0):
            raise ValueError(f"Invalid radial samples for {name}")
        ax.plot(samples[:, 0], samples[:, 1], label=label,
                color=color, linestyle=linestyle, linewidth=2.2)

    if normalization == "equal_cloud_energies":
        reference_energy = metadata["reference_energy"]
        comparison_title = rf"$E_1=E_2\simeq {reference_energy:.2f}\,F_a^2M$"
        caption = ("Each linear cloud is normalized to the saved Omiya reference energy.\n"
                   "Nearest saved nonlinear branch point to the paper's approximate 660; finite-grid comparison.")
    else:
        primary_amplitude, secondary_amplitude = metadata["cloud_amplitudes"]
        comparison_title = "Omiya and linear reference profiles"
        caption = (f"Fixed outer amplitudes: A₁ = {primary_amplitude:.8g}, "
                   f"A₂ = {secondary_amplitude:.8g}\n"
                   "Finite-grid reference comparison; cloud energies are not constrained equal.")
    ax.set_title(
        comparison_title + "\n" +
        rf"$\alpha={metadata['alpha']:g},\quad a/M={metadata['spin']:g}$"
    )
    ax.set_xlabel(r"$r_*/M$", fontsize=14)
    ax.set_ylabel(r"$|R_{lm\omega}|$", fontsize=14)
    ax.set_xlim(-100, 200)
    ax.set_ylim(bottom=0)
    ax.axvline(0, color="0.75", linewidth=0.8, zorder=0)
    ax.grid(alpha=0.16)
    ax.legend(frameon=False, fontsize=11)
    ax.tick_params(direction="in", top=True, right=True)

    # Read the normalization from the data, so the historical fixed-amplitude
    # export can still be plotted with its own accurate caption.
    fig.text(0.5, 0.025, caption, ha="center", fontsize=9, color="0.3")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    args.figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in (".png", ".pdf"):
        destination = args.figure_prefix.with_suffix(extension)
        fig.savefig(destination, dpi=200, bbox_inches="tight")
        print(destination)
    plt.close(fig)


if __name__ == "__main__":
    main()
