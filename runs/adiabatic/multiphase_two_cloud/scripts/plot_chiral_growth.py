"""Plot a complete or partial fixed-BH chiral trajectory and its charge budgets."""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import plotstyle


def plot_chiral_growth(prefix, steps=1000, potential="chiral"):
    trajectory = np.atleast_1d(np.genfromtxt(
        prefix + "-steps%d.csv" % steps, delimiter=",", names=True))
    reference = np.atleast_1d(np.genfromtxt(
        prefix + "-linear-reference.csv", delimiter=",", names=True))
    time = trajectory["time"]
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for index, amplitude in enumerate(("A1", "A2")):
        cloud = reference[reference["cloud"] == index + 1][0]
        initial = cloud["initial_amplitude"]
        axes[0, index].semilogy(time, trajectory[amplitude] / initial, label="Joint " + potential)
        axes[0, index].semilogy(time, np.exp(cloud["linear_growth_rate"] * time),
                              "--", label="Dilute linear reference")
        axes[0, index].set_title(("211", "322")[index])
        axes[0, index].set_ylabel(r"$A(t)/A(0)$")
        axes[1, 0].plot(time, trajectory["log_rate%d" % (index + 1)],
                        label=("211", "322")[index])
    axes[1, 0].set_yscale("symlog", linthresh=1e-12)
    axes[1, 0].set_ylabel(r"$d\ln A/d(t/M)$")
    for charge, error, label in (("energy", "energy_balance_error", "Energy"),
                                 ("angular_momentum", "angular_balance_error", "Angular momentum")):
        scale = np.maximum(abs(trajectory[charge]), abs(trajectory[charge][0]))
        axes[1, 1].plot(time, trajectory[error] / scale, label=label)
    axes[1, 1].set_yscale("symlog", linthresh=1e-12)
    axes[1, 1].set_ylabel("Charge budget error / charge scale")
    for axis in axes.flat:
        axis.set_xlabel(r"$t/M$")
        axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    figure.suptitle("Joint " + potential + " evolution: fixed Kerr background\n"
                   r"$\alpha=0.42$, $\chi=0.99$; exploratory finite-domain trajectory")
    for suffix in (".png", ".pdf"):
        output = prefix + "-growth" + suffix
        figure.savefig(output, dpi=180)
        print(output)
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix")
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--potential", choices=("chiral", "cosine"), default="chiral")
    arguments = parser.parse_args()
    plot_chiral_growth(arguments.prefix, arguments.steps, arguments.potential)
