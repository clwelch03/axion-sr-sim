"""Gather the adaptive Stage E array outputs -> table + Gamma(alpha) plot vs
Leaver, with the 3D DMD point (a=0.9,alpha=0.5 damped, 0.92x) overlaid to show
the 1+1D tool tracks Leaver across the damped branch where 3D could not."""
import numpy as np, glob, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

files = sorted(glob.glob("stage_e_adaptive_[0-9].npz"))
rows = [np.load(f) for f in files]
rows.sort(key=lambda d: float(d["alpha"]))

print(f"{'alpha':>6}{'Leaver G':>13}{'logfit':>9}{'prony':>9}{'freqflt':>9}")
al, Gl, Gff, Glf, Gp = [], [], [], [], []
for d in rows:
    a_ = float(d["alpha"]); G = float(d["G_leaver"])
    g1, g2, g3 = float(d["G_logfit"]), float(d["G_prony"]), float(d["G_freqfilter"])
    print(f"{a_:>6.2f}{G:>13.3e}{g1/G:>8.2f}x{g2/G:>8.2f}x{g3/G:>8.2f}x")
    al.append(a_); Gl.append(G); Glf.append(g1); Gp.append(g2); Gff.append(g3)
al = np.array(al); Gl = np.array(Gl); Gff = np.array(Gff); Glf = np.array(Glf); Gp = np.array(Gp)

fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.plot(al, np.abs(Gl), "k--", marker="s", ms=7, label="Leaver (truth)")
ax.plot(al, np.abs(Gff), "C2", marker="D", ms=6, label="1+1D freqfilter")
ax.plot(al, np.abs(Glf), "C0", marker="o", ms=5, alpha=0.6, label="1+1D logfit")
ax.plot(al, np.abs(Gp), "C1", marker="^", ms=5, alpha=0.6, label="1+1D prony")
ax.scatter([0.50], [5.2050e-04*0.92], color="red", marker="*", s=220,
           label="3D DMD (0.92x, |err|~x2)", zorder=6)
ax.set_yscale("log"); ax.set_xlabel(r"$\alpha = M\mu$"); ax.set_ylabel(r"$|\Gamma|$  (M=1)")
ax.set_title(r"Damped-branch $\Gamma(\alpha)$ at $a=0.9$: 1+1D vs Leaver vs 3D")
ax.legend(); ax.grid(alpha=0.3, which="both")
plt.tight_layout(); plt.savefig("stage_e_adaptive.png", dpi=130)
print("saved stage_e_adaptive.png")
