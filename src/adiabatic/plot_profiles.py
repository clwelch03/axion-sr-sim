from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[2]
data = np.genfromtxt(
    root / "runs/adiabatic/profiles.csv",
    delimiter=",",
    names=True,
    dtype=None,
    encoding="utf-8",
)

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

# Inspect the dominant mode first, keeping both integration branches visible.
for branch, color in (("horizon", "tab:blue"), ("outer", "tab:orange")):
    selected = data[
        (data["branch"] == branch) & (data["l"] == 1) & (data["m"] == 1)
    ]
    selected = np.sort(selected, order="rstar")

    for axis, field in zip(axes, ("R", "Q")):
        axis.plot(
            selected["rstar"], selected[f"{field}_real"],
            color=color, label=f"{branch}: real",
        )
        axis.plot(
            selected["rstar"], selected[f"{field}_imag"],
            color=color, linestyle="--", label=f"{branch}: imaginary",
        )

axes[0].set_ylabel(r"$R_{11}$")
axes[1].set_ylabel(r"$Q_{11}=dR_{11}/dr_*$")
axes[1].set_xlabel(r"$r_*/M$")

for axis in axes:
    axis.axvline(5.0, color="0.6", linewidth=0.8)
    axis.grid(alpha=0.2)

axes[0].legend(ncol=2)
fig.tight_layout()

output_path = root / "plots/adiabatic_profiles.png"
output_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_path, dpi=180)
print(f"Saved {output_path}")