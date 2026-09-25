"""Merge completed cosine derivative shards and plot their comparison to local fits.

Run from any directory using the BHsuper Python environment after all four
tasks finish. The frozen input CSV supplies the expected amplitude coverage.
"""
from pathlib import Path
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from common import plotstyle

CAMPAIGN = Path(__file__).resolve().parents[1]


def main():
    expected = np.atleast_1d(np.genfromtxt(
        CAMPAIGN / "seeds/cosine_branch_20260924.csv", delimiter=",", names=True,
    ))
    parts = [np.atleast_1d(np.genfromtxt(
        CAMPAIGN / f"results/cosine_derivatives_part{index}.csv",
        delimiter=",", names=True,
    )) for index in range(4)]
    data = np.concatenate(parts)
    data = np.sort(data, order="branch_index")
    if not np.array_equal(data["branch_index"], np.arange(1, len(expected) + 1)):
        raise ValueError("missing or duplicate derivative points; wait for all tasks")
    if not np.array_equal(data["amplitude"], expected["amplitude"]):
        raise ValueError("amplitudes differ from the frozen source CSV")
    if not all(np.all(np.isfinite(data[name])) for name in data.dtype.names):
        raise ValueError("nonfinite output; inspect the branch before plotting rates")
    if np.max(data["residual"]) > 1e-10:
        raise ValueError("matching tolerance was not met")

    merged_path = CAMPAIGN / "results/cosine_derivatives.csv"
    np.savetxt(merged_path, np.column_stack([data[name] for name in data.dtype.names]),
               delimiter=",", header=",".join(data.dtype.names), comments="", fmt="%.17g")

    order = np.argsort(data["amplitude"])
    data = data[order]
    amplitudes, energies = data["amplitude"], data["energy"]
    local_derivatives = []
    for index in range(3, len(data) - 3):
        slopes = []
        for width, degree in ((1, 2), (2, 3), (3, 3)):
            selection = slice(index - width, index + width + 1)
            offsets = amplitudes[selection] - amplitudes[index]
            scale = np.max(np.abs(offsets))
            matrix = np.vander(offsets / scale, degree + 1, increasing=True)
            coefficients = np.linalg.lstsq(matrix, energies[selection], rcond=None)[0]
            slopes.append(coefficients[1] / scale)
        local_derivatives.append(slopes)
    local_derivatives = np.array(local_derivatives)
    local_rates = -data["total_flux"][3:-3, None] / local_derivatives
    relative_difference = 100 * (local_derivatives / data["energy_derivative"][3:-3, None] - 1)

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, layout="constrained")
    fig.suptitle("Cosine branch: implicit energy derivatives", fontsize=17, fontweight="bold")
    axes[0].set_title("Fixed spatial grid · preliminary energy and flux numerics", fontsize=11)
    axes[0].plot(amplitudes, data["amplitude_rate"], color="#333333", lw=2,
                 label="Implicit derivative", zorder=5)
    labels = ("3-point quadratic", "5-point cubic", "7-point cubic")
    for column, label in enumerate(labels):
        color = plotstyle.PALETTE[column]
        axes[0].plot(amplitudes[3:-3], local_rates[:, column], linestyle="none",
                     marker=".", ms=4, color=color, label=label)
        axes[1].scatter(amplitudes[3:-3], relative_difference[:, column],
                        s=13, color=color, label=label)
    if np.all(data["amplitude_rate"] > 0):
        axes[0].set_yscale("log")
    axes[0].set_ylabel(r"Amplitude growth $dA_0/d(t/M)$")
    axes[0].legend(frameon=False, fontsize=10)
    axes[1].axhline(0, color="#777777", lw=1)
    axes[1].set_ylabel("Local-fit energy derivative\nrelative to implicit result (%)")
    axes[1].set_xlabel(r"Amplitude $A_0$")
    for axis in axes:
        axis.grid(alpha=0.18)
    fig.supxlabel("Local fits use recomputed energies; dots are individual samples.\n"
                  "Implicit differentiation removes fitting choices, not spatial discretization error.",
                  fontsize=10)
    for extension in ("png", "pdf"):
        fig.savefig(ROOT / f"plots/amplitude_growth_implicit_comparison.{extension}", dpi=180)

    print(f"Merged {len(data)} points: {merged_path}")
    print(f"Amplitude range: {amplitudes[0]:.9g} to {amplitudes[-1]:.9g}")
    print(f"Maximum matching residual: {np.max(data['residual']):.4g}")
    print(f"Energy derivative range: {np.min(data['energy_derivative']):.9g} to {np.max(data['energy_derivative']):.9g}")
    for column, label in enumerate(labels):
        index = np.argmax(np.abs(relative_difference[:, column]))
        print(f"{label}: maximum derivative difference {relative_difference[index, column]:.4f}% at A={amplitudes[3:-3][index]:.9g}")


if __name__ == "__main__":
    main()
