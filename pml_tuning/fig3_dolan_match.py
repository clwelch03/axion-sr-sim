"""Replicate Dolan Fig. 3 with his ACTUAL initial data/parameters (not our
a=0.9 alpha=0.5 validation benchmark): a=0.99, time-symmetric Gaussian
x0=10M, sigma_x=5M, snapshot at t=140M, PML centred around x=-100M,
time series probed at x_1=0. Explicit apples-to-apples comparison to his
Fig. 3, distinct from pml_validation_fig3.png (our own benchmark params).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
import evolve as ev
from radial1d import tortoise_x, horizons

a = 0.99
rp, rm = horizons(a)
# eps_match such that x_match = -100 (Dolan: "a PML centred around x=-100M")
eps_match = brentq(lambda eps: tortoise_x(rp+eps, a) - (-100.0), 1e-10, 1.0)
# probe at x_1 = 0 (his Fig. 3c extraction point)
r_x1 = brentq(lambda r: tortoise_x(r, a) - 0.0, rp+1e-6, 50.0)

common = dict(alpha=1e-6, a=a, N=3500, xmin=-160., xmax=100., T=300.,
              eps_match=eps_match, x0=10., sigx=5., probe_r=r_x1, nsnap=16)
print(f"a={a}  x_match(eps_match={eps_match:.2e})=-100  probe r(x=0)={r_x1:.4f}", flush=True)

print("running no-PML (gmax=0)...", flush=True)
o0 = ev.run(gmax=0.0, **common)
print("running with-PML (gmax=1.0)...", flush=True)
o1 = ev.run(gmax=1.0, **common)

tsnap = 140.
i0 = np.argmin(np.abs(o0["field_t"] - tsnap)); i1 = np.argmin(np.abs(o1["field_t"] - tsnap))
x0g, psi0 = o0["x"], o0["field_psi"][i0]
x1g, psi1 = o1["x"], o1["field_psi"][i1]
xm = o1["x_match"]

fig, axs = plt.subplots(3, 1, figsize=(7, 12))
fig.suptitle("Matching Dolan (2012) Fig. 3 initial data: a=0.99, x0=10M, sigma_x=5M")

axs[0].plot(x0g, psi0.real, color="r", label="(i) without PML")
axs[0].plot(x1g, psi1.real, "b--", label="(ii) with PML")
axs[0].axvline(xm, color="m", ls=":", label="x_match")
axs[0].set_xlabel("x/M"); axs[0].set_ylabel("Re(psi)")
axs[0].set_title(f"(a) snapshot at t={o1['field_t'][i1]:.0f}M"); axs[0].legend()
axs[0].set_xlim(-150, 100); axs[0].set_ylim(-0.1, 0.06)

eps = 1e-300
axs[1].semilogy(x0g, np.abs(psi0.real) + eps, color="r", label="(i) without PML")
axs[1].semilogy(x1g, np.abs(psi1.real) + eps, "b--", label="(ii) with PML")
axs[1].axvline(xm, color="m", ls=":")
axs[1].set_xlabel("x/M"); axs[1].set_ylabel("|Re(psi)|")
axs[1].set_title("(b) same, log scale"); axs[1].legend()
axs[1].set_xlim(-150, 100); axs[1].set_ylim(1e-12, 1)

axs[2].semilogy(o0["t"], np.abs(o0["psi"].real) + eps, color="r", label="(i) without PML")
axs[2].semilogy(o1["t"], np.abs(o1["psi"].real) + eps, "b--", label="(ii) with PML")
axs[2].set_xlabel("t/M"); axs[2].set_ylabel("|Re(psi(x=0))|")
axs[2].set_title("(c) probe time series at x_1=0"); axs[2].legend()
axs[2].set_xlim(0, 300); axs[2].set_ylim(1e-10, 1)

plt.tight_layout()
plt.savefig("pml_validation_fig3_DOLAN_MATCH.png", dpi=130)
print("saved pml_validation_fig3_DOLAN_MATCH.png", flush=True)
