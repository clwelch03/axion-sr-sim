"""Replicate Dolan Fig. 3: validate the PML by comparing no-PML (gmax=0) vs
with-PML (gmax=<best>) evolutions of the *same* initial data.
  (a) snapshot Re(psi(x)) at fixed t, linear scale
  (b) same snapshot, log scale
  (c) probe time series (outside the PML) vs t, log scale, long enough to see
      the no-PML run's wall-reflection arrive and split off from with-PML.
Uses the current best PML config (eps_match=0.002, xmin=-120) from scan_pml.py.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
import evolve as ev

alpha, a = 0.50, 0.90
common = dict(alpha=alpha, a=a, N=4000, xmin=-120., xmax=220., T=600.,
              eps_match=0.002, nsnap=13)

print("running no-PML (gmax=0)...", flush=True)
o0 = ev.run(gmax=0.0, **common)
print("running with-PML (gmax=1.5)...", flush=True)
o1 = ev.run(gmax=1.5, **common)

tsnap = 150.
i0 = np.argmin(np.abs(o0["field_t"] - tsnap)); i1 = np.argmin(np.abs(o1["field_t"] - tsnap))
x0, psi0 = o0["x"], o0["field_psi"][i0]
x1, psi1 = o1["x"], o1["field_psi"][i1]
xm = o1["x_match"]

fig, axs = plt.subplots(3, 1, figsize=(7, 12))

axs[0].plot(x0, psi0.real, color="r", label="(i) without PML")
axs[0].plot(x1, psi1.real, "b--", label="(ii) with PML")
axs[0].axvline(xm, color="m", ls=":", label="x_match")
axs[0].set_xlabel("x/M"); axs[0].set_ylabel("Re(psi)")
axs[0].set_title(f"(a) snapshot at t={o1['field_t'][i1]:.0f}M"); axs[0].legend()

eps = 1e-300
axs[1].semilogy(x0, np.abs(psi0.real) + eps, color="r", label="(i) without PML")
axs[1].semilogy(x1, np.abs(psi1.real) + eps, "b--", label="(ii) with PML")
axs[1].axvline(xm, color="m", ls=":")
axs[1].set_xlabel("x/M"); axs[1].set_ylabel("|Re(psi)|")
axs[1].set_title("(b) same, log scale"); axs[1].legend()

axs[2].semilogy(o0["t"], np.abs(o0["psi"].real) + eps, color="r", label="(i) without PML")
axs[2].semilogy(o1["t"], np.abs(o1["psi"].real) + eps, "b--", label="(ii) with PML")
axs[2].set_xlabel("t/M"); axs[2].set_ylabel(f"|Re(psi(r={o1['probe_r']:.1f}))|")
axs[2].set_title("(c) probe time series (outside PML)"); axs[2].legend()

plt.tight_layout()
plt.savefig("pml_validation_fig3.png", dpi=130)
print("saved pml_validation_fig3.png", flush=True)
