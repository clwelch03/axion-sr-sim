"""Stage D floor test (cheap decisive test before committing to a 1e5 M run).

Attempt 1 failed (Gamma ~ -6e-6, wrong SIGN, 40x too big) for two reasons:
PML was only 0.83*lambda wide AND the run was 5-10x too short vs Dolan's
t_max=1e5 M. This test isolates the FIRST: with a properly-sized PML
(2.82*lambda, following our own damped-scan calibration width >= 2.5*lambda,
lambda=2pi/|E-m*OmH|=251 M here), how low is the systematic Gamma floor?

Params match Dolan Fig.13 where possible: dx = M/6 (=0.16667), a=0.99,
alpha=0.42 (the growth peak, Leaver Gamma=+1.52e-7). Moderate T=1e4 (not the
full 1e5) -- enough to READ the floor magnitude cheaply.

Decision rule:
  - floor drops to <~1e-6 (from attempt-1's ~1e-5)  -> good PML works, a
    long (~1e5 M) run can plausibly resolve +1.5e-7; proceed to full run.
  - floor stuck ~1e-5                                -> pure-1D-numpy can't
    reach 1.5e-7 at the peak; need frequency-filter / C port / longer still.
"""
import numpy as np, sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
import evolve as ev

alpha, a = 0.42, 0.99
w = complex(leaver.solve(alpha, a))
print(f"Leaver target: E={w.real:.6f}  Gamma={w.imag:+.4e}  (want SIGN=+)", flush=True)

# dx = M/6 = 0.16667 exactly (Dolan Fig.13); xmin=-750 -> width 2.82*lambda
o = ev.run(alpha=alpha, a=a, N=6000, xmin=-750., xmax=250., T=10000.,
           eps_match=0.002, gmax=1.5, nrec=8000)
print(f"probe r={o['probe_r']:.2f}  dx={o['dx']:.4f}  x_match={o['x_match']:.2f}"
      f"  PML width={o['x_match']-(-750):.0f} = {(o['x_match']-(-750))/251.3:.2f} lambda", flush=True)

t, psi = o["t"], o["psi"]
E1, G1 = ev.logfit(t, psi)
print(f"logfit : E={E1:.6f}  Gamma={G1:+.4e}  ({G1/w.imag:+.2f}x)", flush=True)
for tmf in [0.1, 0.35, 0.6]:
    try:
        E2, G2, om, c = ev.prony(t, psi, p=12, tmin_frac=tmf, anchor=w, full=True)
        print(f"prony (tmin_frac={tmf}): E={E2:.6f} Gamma={G2:+.4e}  ({G2/w.imag:+.2f}x)"
              f"  |nearest-Leaver Gamma|={abs(G2):.2e}", flush=True)
    except Exception as e:
        print(f"prony (tmin_frac={tmf}) failed: {e}", flush=True)
np.savez("stage_d_floortest.npz", t=t, psi=psi, w_leaver=w)
print("saved stage_d_floortest.npz", flush=True)
