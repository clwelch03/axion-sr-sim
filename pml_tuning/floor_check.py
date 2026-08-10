"""Check for an extraction-window floor: plot log|psi(t)| (cloud probe) and
log|F_E(t)| (flux at x_match) over the FULL T=3000 production window at the
current best PML config, to see if the late-time decay slope flattens
(a floor would explain the logfit/flux agreement at ~0.95x Leaver despite
further PML widening not helping -- Sec. 5c/8 of README.md).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
import evolve as ev

alpha, a = 0.50, 0.90
o = ev.run(alpha=alpha, a=a, N=4000, xmin=-120., xmax=220., T=3000.,
           eps_match=0.002, gmax=1.5)
w = o["w_ref"]
t = o["t"]; psi = o["psi"]
Gflux, F_E = ev.flux_gamma(t, o["phi_flux"], o["phidot_flux"], o["phix_flux"], o["flux_r"], a)

# local slope of log|signal| in sliding windows, to see if it flattens late
def local_slope(t, y, frac=0.05):
    n = len(t); w_ = max(int(n*frac), 10)
    ts, sl = [], []
    for i in range(0, n-w_, w_//2):
        tt = t[i:i+w_]; yy = np.log(np.abs(y[i:i+w_])+1e-300)
        s = np.polyfit(tt, yy, 1)[0]
        ts.append(tt.mean()); sl.append(s)
    return np.array(ts), np.array(sl)

ts_psi, sl_psi = local_slope(t, psi)
ts_F, sl_F = local_slope(t, F_E)  # this slope is ~2*Gamma

fig, axs = plt.subplots(3, 1, figsize=(8, 12))

axs[0].semilogy(t, np.abs(psi)+1e-300, label="|psi| (cloud probe)")
axs[0].axvline(0.4*t[-1], color="gray", ls=":", label="logfit tmin_frac=0.4")
axs[0].set_xlabel("t/M"); axs[0].set_ylabel("|psi|"); axs[0].legend()
axs[0].set_title("full-window amplitude decay")

axs[1].plot(ts_psi, sl_psi, "o-", label="local slope of log|psi| (= Gamma)")
axs[1].axhline(w.imag, color="k", ls="--", label=f"Leaver Gamma={w.imag:.3e}")
axs[1].axvline(0.4*t[-1], color="gray", ls=":")
axs[1].set_xlabel("t/M"); axs[1].set_ylabel("local d ln|psi|/dt")
axs[1].set_title("sliding-window slope: flattening = floor"); axs[1].legend()

axs[2].plot(ts_F, sl_F/2.0, "o-", color="C1", label="local slope of log|F_E|/2 (= Gamma_flux)")
axs[2].axhline(w.imag, color="k", ls="--", label=f"Leaver Gamma={w.imag:.3e}")
axs[2].axvline(0.4*t[-1], color="gray", ls=":")
axs[2].set_xlabel("t/M"); axs[2].set_ylabel("local (d ln|F_E|/dt)/2")
axs[2].set_title("same, from the horizon-interface flux"); axs[2].legend()

plt.tight_layout()
plt.savefig("floor_check.png", dpi=130)
print(f"Leaver Gamma = {w.imag:+.4e}")
print(f"logfit (t>0.4T) Gamma = {ev.logfit(t,psi)[1]:+.4e}")
print(f"flux   (t>0.4T) Gamma = {Gflux:+.4e}")
print("last 5 sliding-window psi slopes:", sl_psi[-5:])
print("last 5 sliding-window flux slopes/2:", (sl_F/2.0)[-5:])
print("saved floor_check.png")
