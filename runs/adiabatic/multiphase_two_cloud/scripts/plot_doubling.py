"""Compare the joint trajectory with linear growth and an optional Omiya trajectory."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import plotstyle


def plot_doubling(prefix, steps=64, spin=None, alpha=None, omiya_path=None):
    """TODO: describe normalized amplitudes and fractional reference differences."""
    trajectory = np.atleast_1d(np.genfromtxt(
        str(prefix) + "-steps%d.csv" % steps, delimiter=",", names=True))
    reference = np.atleast_1d(np.genfromtxt(
        str(prefix) + "-linear-reference.csv", delimiter=",", names=True))
    time = trajectory["time"]
    omiya = None
    if omiya_path is not None:
        omiya = np.atleast_1d(np.genfromtxt(omiya_path, delimiter=",", names=True))
        if len(omiya) != len(time) or not np.allclose(omiya["time"], time, rtol=1e-12, atol=0):
            raise ValueError("Omiya trajectory must use the joint trajectory's saved times")
        for amplitude in ("A1", "A2"):
            if not np.all(np.isfinite(omiya[amplitude])) or np.any(omiya[amplitude] <= 0):
                raise ValueError("Omiya amplitudes must be finite and positive")
            if not np.isclose(omiya[amplitude][0], trajectory[amplitude][0], rtol=1e-12, atol=0):
                raise ValueError("Omiya and joint initial amplitudes must match")
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    comparison = [time]
    columns = ["time"]
    for index, amplitude in enumerate(("A1", "A2")):
        cloud = reference[reference["cloud"] == index + 1][0]
        initial_amplitude = cloud["initial_amplitude"]
        # An amplitude grows with Im(omega); the energy exponent is twice that.
        linear = initial_amplitude * np.exp(cloud["linear_growth_rate"] * time)
        frozen = initial_amplitude * np.exp(cloud["initial_log_rate"] * time)
        evolved = trajectory[amplitude]
        axes[0, index].plot(time, evolved / initial_amplitude, label="Joint cosine evolution")
        axes[0, index].plot(time, linear / initial_amplitude, "--", label="Linear eigenmode")
        if omiya is None:
            axes[0, index].plot(time, frozen / initial_amplitude, ":", label="Initial E/J growth rate")
            axes[1, index].plot(time, evolved / linear, label="Joint / linear eigenmode")
            axes[1, index].plot(time, evolved / frozen, label="Joint / initial E/J rate")
        else:
            reference_amplitude = omiya[amplitude]
            axes[0, index].plot(time, reference_amplitude / initial_amplitude, ":",
                               label="Omiya approximation")
            axes[1, index].plot(time, evolved / linear, color="C0",
                               label="Joint / linear")
            axes[1, index].plot(time, reference_amplitude / linear, color="C2",
                               linestyle=":", label="Omiya / linear")
        axes[0, index].set_title(("Primary: 211", "Secondary: 322")[index])
        axes[0, index].set_ylabel(r"$A(t)/A(0)$")
        axes[0, index].ticklabel_format(axis="y", style="plain", useOffset=False)
        axes[1, index].set_ylabel(r"$A(t)/A_{\rm linear}(t)$" if omiya is not None
                                 else "Amplitude / reference amplitude")
        axes[1, index].ticklabel_format(axis="y", style="plain", useOffset=False)
        comparison.extend((evolved, linear, frozen, evolved / linear - 1, evolved / frozen - 1))
        columns.extend(amplitude + suffix for suffix in (
            "", "_linear", "_initial_rate", "_relative_to_linear", "_relative_to_initial_rate"))
        if omiya is not None:
            comparison.extend((reference_amplitude, reference_amplitude / linear - 1,
                               evolved / reference_amplitude - 1))
            columns.extend(amplitude + suffix for suffix in (
                "_omiya", "_omiya_relative_to_linear", "_relative_to_omiya"))
    for axis in axes.flat:
        axis.set_xlabel(r"$t/M$")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
        axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    title = "Primary-amplitude doubling on a fixed Kerr background\n"
    initial_conditions = []
    if spin is not None:
        initial_conditions.append(r"$a_*(0)=J_{\rm BH}/M^2=%.4g$" % spin)
    if alpha is not None:
        initial_conditions.append(r"$\alpha(0)=M\mu=%.4g$" % alpha)
    if initial_conditions:
        title += "Initial conditions: " + ",  ".join(initial_conditions) + "\n"
    figure.suptitle(title + "Exploratory finite-grid evolution")
    for suffix in (".png", ".pdf"):
        output = str(prefix) + "-comparison" + suffix
        figure.savefig(output, dpi=180)
        print(output)
    plt.close(figure)
    np.savetxt(str(prefix) + "-comparison.csv", np.column_stack(comparison),
               delimiter=",", header=",".join(columns), comments="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--spin", type=float, help="Initial dimensionless BH spin J_BH/M^2")
    parser.add_argument("--alpha", type=float, help="Initial gravitational coupling M*mu")
    parser.add_argument("--omiya", type=Path, help="Matched Omiya-approximation trajectory CSV")
    arguments = parser.parse_args()
    plot_doubling(arguments.prefix, arguments.steps, arguments.spin, arguments.alpha,
                  arguments.omiya)
