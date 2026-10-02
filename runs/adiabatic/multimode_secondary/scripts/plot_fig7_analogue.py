"""Fig. 7 analogue: nonlinear and linear |R_11|, |R_22| at E1 = E2.

PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig7_analogue.py alpha030 660

Reads results/<tag>_cosine_branch.csv, results/<tag>_cosine_secondary_energy.csv,
and work/<tag>_profiles.csv (export_profiles.jl <tag> <energy>). Writes
plots/multimode_<tag>_fig7_analogue.png. Linear profiles are the smallest saved
amplitude, rescaled to the same energy. Omiya et al. 2211.01949 Fig. 7 is at
Mμ = 0.3, a/M = 0.99, E1 = E2 ~ 660 Fa² M.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import plotstyle

tag, target_energy = sys.argv[1], float(sys.argv[2])
campaign = Path(__file__).resolve().parents[1]
root = campaign.parents[2]
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.size": 11})
primary_color, secondary_color = plotstyle.PALETTE[0], plotstyle.PALETTE[1]

primary = np.genfromtxt(campaign / f"results/{tag}_cosine_branch.csv", delimiter=",", names=True)
secondary = np.genfromtxt(campaign / f"results/{tag}_cosine_secondary_energy.csv",
                          delimiter=",", names=True)
profiles = np.genfromtxt(campaign / f"work/{tag}_profiles.csv", delimiter=",", names=True)
if not np.allclose(primary["amplitude"][:len(secondary)], secondary["amplitude"]):
    raise ValueError("primary and secondary tables are not aligned")

target_index = np.argmin(np.abs(primary["energy"][:len(secondary)] - target_energy))
target_amplitude = primary["amplitude"][target_index]
energy = primary["energy"][target_index]


def profile_at(value):
    rows = profiles[np.isclose(profiles["amplitude"], value, rtol=1e-8)]
    if len(rows) == 0:
        raise ValueError(f"no exported profile at A1 = {value}")
    return rows[np.argsort(rows["tortoise_radius"])]


nonlinear, linear = profile_at(target_amplitude), profile_at(primary["amplitude"][0])
# Primary R is physical; secondary R is per unit outer amplitude.
linear_primary_scale = np.sqrt(energy / primary["energy"][0])
nonlinear_secondary_scale = np.sqrt(energy / secondary["energy"][target_index])
linear_secondary_scale = np.sqrt(energy / secondary["energy"][0])

fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
axis.plot(nonlinear["tortoise_radius"], nonlinear["primary_abs_R11"],
          color=primary_color, lw=2, label=r"nonlinear $l=m=1$")
axis.plot(linear["tortoise_radius"], linear_primary_scale * linear["primary_abs_R11"],
          color=primary_color, lw=1.5, ls="--", label=r"linear $l=m=1$")
axis.plot(nonlinear["tortoise_radius"], nonlinear_secondary_scale * nonlinear["secondary_abs_R22"],
          color=secondary_color, lw=2, label=r"nonlinear $l=m=2$")
axis.plot(linear["tortoise_radius"], linear_secondary_scale * linear["secondary_abs_R22"],
          color=secondary_color, lw=1.5, ls="--", label=r"linear $l=m=2$")
axis.set_xlim(-100, 200)
axis.set_xlabel(r"tortoise radius $r_*/M$")
axis.set_ylabel(r"$|R_{lm}|$")
axis.legend(loc="upper right", frameon=False)
axis.grid(True, which="major", color="0.92", lw=0.6)
axis.set_title(rf"$E_1=E_2\approx{energy:.0f}\,F_a^2M$ ($A_1={target_amplitude:.3g}$), {tag}",
               fontsize=11)
output = root / f"plots/multimode_{tag}_fig7_analogue.png"
fig.savefig(output, dpi=160)

# Peak values for comparison with the paper's figure.
for name, column, scale_nonlinear, scale_linear in (
        ("l=m=1", "primary_abs_R11", 1.0, linear_primary_scale),
        ("l=m=2", "secondary_abs_R22", nonlinear_secondary_scale, linear_secondary_scale)):
    peak_nonlinear = np.argmax(nonlinear[column])
    peak_linear = np.argmax(linear[column])
    print(f"{name}: nonlinear peak {scale_nonlinear * nonlinear[column][peak_nonlinear]:.4f} "
          f"at r*={nonlinear['tortoise_radius'][peak_nonlinear]:.1f}; linear peak "
          f"{scale_linear * linear[column][peak_linear]:.4f} at r*={linear['tortoise_radius'][peak_linear]:.1f}")
print("wrote", output)
