"""Fig. 6 analogue: nonperturbative evolution of the two cloud energies.

PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig6.py

Reads results/two_cloud_evolution_alpha042.csv and _alpha030.csv (written by
src/adiabatic/drivers/two_cloud_evolution.jl). Writes
plots/multimode_fig6_analogue.png. Omiya et al. 2211.01949 Fig. 6: a/M = 0.99,
Mμ = 0.42 and 0.3, E1(0) ~ 0.27 and 12 Fa² M, E2(0) = 1e-8 Fa² M, time in units
of the linear growth rate ω_{1,I}. Only the nonperturbative curves are drawn.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import plotstyle

campaign = Path(__file__).resolve().parents[1]
root = campaign.parents[2]
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.size": 11})
primary_color, secondary_color = plotstyle.PALETTE[0], plotstyle.PALETTE[1]

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
for axis, (label, alpha, top) in zip(axes, (("alpha042", 0.42, 2000), ("alpha030", 0.3, 4000))):
    data = np.genfromtxt(campaign / f"results/two_cloud_evolution_{label}.csv",
                         delimiter=",", names=True)
    axis.semilogy(data["time"], data["primary_energy"], color=primary_color, lw=2,
                  label=r"nonperturbative $l=m=1$")
    axis.semilogy(data["time"], data["secondary_energy"], color=secondary_color, lw=2, ls=":",
                  label=r"nonperturbative $l=m=2$")
    axis.set_xlim(0, 10)
    axis.set_ylim(1, top)
    axis.set_xlabel(r"$\omega_{1,I}\,t$")
    axis.set_ylabel(r"$E/(F_a^2M)$")
    axis.set_title(rf"$a/M=0.99,\ \mu M={alpha}$")
    axis.grid(True, which="major", color="0.92", lw=0.6)
    axis.legend(loc="lower right", frameon=False)

output = root / "plots/multimode_fig6_analogue.png"
fig.savefig(output, dpi=160)
print("wrote", output)
