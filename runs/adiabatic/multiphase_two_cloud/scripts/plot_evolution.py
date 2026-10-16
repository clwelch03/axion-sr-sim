"""TODO: describe the finite-domain pilot plots and comparison inputs."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import plotstyle  # Apply the repository's plotting palette.


def plot_evolution(prefix, output):
    """TODO: describe amplitude changes, charges, and boundary-loss budget errors."""
    fine = np.genfromtxt(str(prefix) + "-steps8.csv", delimiter=",", names=True)
    coarse = np.genfromtxt(str(prefix) + "-steps4.csv", delimiter=",", names=True)
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for table, marker, label in ((coarse, "o", "4 steps"), (fine, ".", "8 steps")):
        time = table["time"]
        for index, amplitude in enumerate(("A1", "A2")):
            axes[0, index].plot(time, 100 * (table[amplitude] / table[amplitude][0] - 1),
                                marker=marker, label=label)
        for index, (charge, error) in enumerate((("energy", "energy_balance_error"),
                                                ("angular_momentum", "angular_balance_error"))):
            axes[1, index].plot(time, table[error] / table[charge][0],
                                marker=marker, label=label)
    axes[0, 0].set_title("First cloud: 211")
    axes[0, 1].set_title("Second cloud: 322")
    axes[0, 0].set_ylabel("Amplitude change (%)")
    axes[0, 1].set_ylabel("Amplitude change (%)")
    axes[1, 0].set_ylabel(r"$(E-E_0+\int F_E\,dt)/E_0$")
    axes[1, 1].set_ylabel(r"$(J-J_0+\int F_J\,dt)/J_0$")
    for axis in axes.flat:
        axis.set_xlabel(r"$t/M$")
        axis.grid(alpha=0.2)
        axis.legend()
        axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    figure.suptitle("Two-cloud cosine evolution on a fixed Kerr background\n"
                   "Finite-domain pilot; spatial convergence not established")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    plot_evolution(arguments.prefix, arguments.output)
