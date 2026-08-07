"""Adaptive Stage E: Gamma(alpha) sweep at a=0.9 vs Leaver, ONE alpha per
SLURM array task. The fixed-(T, window) sweep only worked at alpha=0.5 because
|Gamma| spans ~6 orders across the band (superradiant +1.5e-8 at 0.3 up to
-8e-3 at 0.6); a single run length can't serve all. Here every knob is derived
per-alpha from the Leaver rate:

  lambda = 2pi/|E - m*OmH|        (near-horizon wavelength)
  xmin   = min(-120, x_match - 2.7*lambda)   (PML width >= 2.5 lambda rule)
  T      = clip(5/|Gamma|, 2500, 6e4)        (span a few e-folds, capped for cost)
  dx     = 0.09 fixed  ->  N = (xmax-xmin)/dx

Primary estimator = freqfilter (Dolan's, validated 1.000x at alpha=0.5): it
reconstructs the QBS mode and log-regresses the CENTRAL 50% of the run, which
the adaptive T places in the good zone (mode established, not yet decayed to
the ~1e-8 floor). logfit/prony reported as secondary.
"""
import numpy as np, sys, os
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
import evolve as ev
from radial1d import horizons, tortoise_x

ALPHAS = [0.40, 0.45, 0.50, 0.55, 0.60]
a = 0.90
idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 2))
alpha = ALPHAS[idx]

w = complex(leaver.solve(alpha, a)); E, G = w.real, w.imag
rp, rm = horizons(a); OmH = a/(2*rp)
lam = 2*np.pi/abs(E - OmH)
xm = tortoise_x(rp + 0.002, a)
xmin = min(-120.0, xm - 2.7*lam)
T = float(np.clip(5.0/abs(G), 2500.0, 60000.0))
xmax = 220.0; dx = 0.09
N = int((xmax - xmin)/dx)
nrec = int(min(16000, max(6000, T/5)))
print(f"[task {idx}] alpha={alpha} a={a}  Leaver E={E:.6f} Gamma={G:+.4e}", flush=True)
print(f"  lambda={lam:.1f}  xmin={xmin:.0f} (PML {(xm-xmin)/lam:.2f} lambda)"
      f"  T={T:.0f}  N={N}  nrec={nrec}", flush=True)

o = ev.run(alpha=alpha, a=a, N=N, xmin=xmin, xmax=xmax, T=T,
           eps_match=0.002, gmax=1.5, nrec=nrec)
t, psi = o["t"], o["psi"]

E1, G1 = ev.logfit(t, psi)
try:
    E2, G2 = ev.prony(t, psi, anchor=w)
except Exception as e:
    E2, G2 = np.nan, np.nan; print("  prony failed:", e, flush=True)
E3, G3 = ev.freqfilter(t, psi, E_peak=E, kbins=5.0)

print(f"  Leaver     Gamma={G:+.4e}", flush=True)
print(f"  logfit     Gamma={G1:+.4e}  ({G1/G:+.3f}x)  E={E1:.6f}", flush=True)
print(f"  prony      Gamma={G2:+.4e}  ({G2/G:+.3f}x)  E={E2:.6f}", flush=True)
print(f"  freqfilter Gamma={G3:+.4e}  ({G3/G:+.3f}x)  E={E3:.6f}", flush=True)

np.savez(f"stage_e_adaptive_{idx}.npz",
         alpha=alpha, E_leaver=E, G_leaver=G,
         G_logfit=G1, G_prony=G2, G_freqfilter=G3,
         E_logfit=E1, E_prony=E2, E_freqfilter=E3,
         lam=lam, xmin=xmin, T=T, N=N)
print(f"  saved stage_e_adaptive_{idx}.npz", flush=True)
