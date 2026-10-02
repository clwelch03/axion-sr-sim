"""Plot the averaged secondary branch and its environment.

PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_secondary.py

Reads results/cosine_secondary.csv and work/profiles.csv (from export_profiles.jl).
Writes plots/multimode_secondary_frequency.png and
plots/multimode_secondary_environment.png.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import plotstyle

campaign = Path(__file__).resolve().parents[1]
root = campaign.parents[2]
plots = root / "plots"
plots.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.size": 11})

primary_color, secondary_color = plotstyle.PALETTE[0], plotstyle.PALETTE[1]
guide_style = dict(color="0.35", lw=1.0, ls="--")
# Above this amplitude the saved primary is not step-converged (campaign README).
unconverged_amplitude = 9.0

branch = np.genfromtxt(campaign / "results/cosine_secondary.csv", delimiter=",", names=True)
linear_real, linear_imag = 0.415628403633, 6.370261e-11  # A1 = 0 solve (README)
amplitude = branch["amplitude"]

# Figure 1: omega2 along the branch.
fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.2), sharex=True, constrained_layout=True)

axes[0].loglog(amplitude, linear_real - branch["omega2_real"], color=secondary_color, lw=2)
axes[0].loglog(amplitude, 0.0101 * amplitude**2, **guide_style)
axes[0].text(0.003, 1.5e-6, r"$0.0101\,A_1^2$", color="0.35")
axes[0].set_ylabel(r"$-\,\delta\,\mathrm{Re}\,M\omega_2$")

axes[1].loglog(amplitude, branch["omega2_imag"] / linear_imag, color=secondary_color, lw=2)
axes[1].loglog(amplitude, 1 + 42 * amplitude**2, **guide_style)
axes[1].text(0.25, 1.2, r"$1 + 42\,A_1^2$", color="0.35")
axes[1].set_ylabel(r"$\Gamma_2 / \Gamma_2^{\mathrm{linear}}$")

growing_primary = branch["omega1_imag"] > 0
axes[2].loglog(amplitude, branch["omega2_imag"], color=secondary_color, lw=2,
               label=r"secondary $\Gamma_2$")
axes[2].loglog(amplitude[growing_primary], branch["omega1_imag"][growing_primary],
               color=primary_color, lw=2, label=r"primary $\Gamma_1$ (where $>0$)")
axes[2].text(0.002, 9e-11, r"$\Gamma_2$", color="0.15")
axes[2].text(0.002, 2.1e-7, r"$\Gamma_1$", color="0.15")
axes[2].legend(loc="lower right", frameon=False)
axes[2].set_ylabel(r"$M\,\Gamma$")
axes[2].set_xlabel(r"primary outer amplitude $A_1$")

for axis in axes:
    axis.axvspan(unconverged_amplitude, amplitude.max(), color="0.9", lw=0)
    axis.grid(True, which="major", color="0.92", lw=0.6)
axes[1].text(unconverged_amplitude * 0.95, 3.0, "primary not\nstep-converged",
             ha="right", va="center", color="0.4", fontsize=9)
axes[0].set_title(r"Averaged secondary on the cosine primary, $a=0.99$, $M\mu=0.42$",
                  fontsize=11)
fig.savefig(plots / "multimode_secondary_frequency.png", dpi=160)

# Figure 2: the secondary's environment.
profiles = np.genfromtxt(campaign / "work/profiles.csv", delimiter=",", names=True)
selected = np.unique(profiles["amplitude"])
shades = plt.cm.Blues(np.linspace(0.35, 0.95, len(selected)))

fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True, constrained_layout=True)
for shade, value in zip(shades, selected):
    rows = profiles[profiles["amplitude"] == value]
    rows = rows[np.argsort(rows["tortoise_radius"])]
    label = rf"$A_1={value:.3g}$"
    axes[0].semilogy(rows["tortoise_radius"], rows["mass_deficit_equator"],
                     color=shade, lw=1.8, label=label)
    axes[1].plot(rows["tortoise_radius"],
                 rows["secondary_abs_R22"] / rows["secondary_abs_R22"].max(),
                 color=shade, lw=1.8, label=label)

axes[0].axhline(1.0, **guide_style)
axes[0].text(148, 1.15, r"$\langle\cos\phi_1\rangle = 0$", ha="right", color="0.35")
axes[0].set_ylim(1e-7, 3)
axes[0].set_ylabel(r"$\langle 1-\cos\phi_1\rangle$ on the equator")
axes[0].legend(loc="lower left", frameon=False, fontsize=9, ncol=2)
axes[1].set_ylabel(r"$|R_{22}|$, normalized to its peak")
axes[1].set_xlabel(r"tortoise radius $r_*/M$")
axes[1].set_xlim(-30, 150)
for axis in axes:
    axis.grid(True, which="major", color="0.92", lw=0.6)
axes[0].set_title("Primary's mass deficit and the secondary's radial profile", fontsize=11)
fig.savefig(plots / "multimode_secondary_environment.png", dpi=160)
print("wrote", plots / "multimode_secondary_frequency.png")
print("wrote", plots / "multimode_secondary_environment.png")

# Figure 3: growth rates against the primary's energy (the paper's axis).
primary_branch = np.genfromtxt(root / "runs/adiabatic/fig2_20260922/results/cosine_branch.csv",
                               delimiter=",", names=True)
secondary_energy = np.genfromtxt(campaign / "results/cosine_secondary_energy.csv",
                                 delimiter=",", names=True)
if not np.allclose(primary_branch["amplitude"], branch["amplitude"]):
    raise ValueError("primary and secondary tables are not aligned")
primary_energy = primary_branch["energy"]
unconverged_energy = np.interp(unconverged_amplitude, amplitude, primary_energy)

fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
axis.loglog(primary_energy, branch["omega2_imag"], color=secondary_color, lw=2,
            label=r"secondary $\Gamma_2$")
axis.loglog(primary_energy[growing_primary], branch["omega1_imag"][growing_primary],
            color=primary_color, lw=2, label=r"primary $\Gamma_1$ (where $>0$)")
axis.axvspan(unconverged_energy, primary_energy.max(), color="0.9", lw=0)
axis.text(unconverged_energy * 0.95, 1e-10, "primary not\nstep-converged",
          ha="right", va="bottom", color="0.4", fontsize=9)
axis.set_xlabel(r"primary energy $E_1/(F_a^2 M)$")
axis.set_ylabel(r"$M\,\Gamma$")
axis.legend(loc="upper left", frameon=False)
axis.grid(True, which="major", color="0.92", lw=0.6)
axis.set_title(r"Growth rates against primary energy, $a=0.99$, $M\mu=0.42$", fontsize=11)
fig.savefig(plots / "multimode_secondary_growth_vs_energy.png", dpi=160)

# Figure 4: Fig. 7 analogue. Radial profiles at E1 = E2 near 1000 Fa² M,
# nonlinear (solid) against linear profiles scaled to the same energy (dashed).
target_index = np.argmin(np.abs(primary_energy - 1000))
target_amplitude, target_energy = amplitude[target_index], primary_energy[target_index]
linear_amplitude = amplitude[0]   # smallest saved amplitude, linear to ~1e-6


def profile_at(value):
    rows = profiles[np.isclose(profiles["amplitude"], value, rtol=1e-8)]
    return rows[np.argsort(rows["tortoise_radius"])]


nonlinear, linear = profile_at(target_amplitude), profile_at(linear_amplitude)
# The primary's R is physical; rescale the linear one to the same energy.
linear_primary_scale = np.sqrt(target_energy / primary_energy[0])
# The secondary's R is per unit outer amplitude; scale each to E2 = target_energy.
nonlinear_secondary_scale = np.sqrt(target_energy / secondary_energy["energy"][target_index])
linear_secondary_scale = np.sqrt(target_energy / secondary_energy["energy"][0])

fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
axis.plot(nonlinear["tortoise_radius"], nonlinear["primary_abs_R11"],
          color=primary_color, lw=2, label=r"nonlinear $l=m=1$")
axis.plot(linear["tortoise_radius"], linear_primary_scale * linear["primary_abs_R11"],
          color=primary_color, lw=1.5, ls="--", label=r"linear $l=m=1$")
axis.plot(nonlinear["tortoise_radius"], nonlinear_secondary_scale * nonlinear["secondary_abs_R22"],
          color=secondary_color, lw=2, label=r"nonlinear $l=m=2$")
axis.plot(linear["tortoise_radius"], linear_secondary_scale * linear["secondary_abs_R22"],
          color=secondary_color, lw=1.5, ls="--", label=r"linear $l=m=2$")
axis.set_xlim(-30, 150)
axis.set_xlabel(r"tortoise radius $r_*/M$")
axis.set_ylabel(r"$|R_{lm}|$")
axis.legend(loc="upper right", frameon=False)
axis.grid(True, which="major", color="0.92", lw=0.6)
axis.set_title(rf"$E_1=E_2\approx{target_energy:.0f}\,F_a^2M$ ($A_1={target_amplitude:.3g}$), "
               r"$a=0.99$, $M\mu=0.42$", fontsize=11)
fig.savefig(plots / "multimode_secondary_fig7_analogue.png", dpi=160)
print("wrote", plots / "multimode_secondary_growth_vs_energy.png")
print("wrote", plots / "multimode_secondary_fig7_analogue.png")
