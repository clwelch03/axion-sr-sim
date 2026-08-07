"""Stage D DEFINITIVE run: the actual deliverable -- a scalar superradiant
GROWTH rate from the 1+1D time-domain code, at the a=0.99 growth peak
(alpha=0.42, Leaver Gamma=+1.5183e-07).

All parameters now justified (not guessed):
  - dx = M/6 = 0.16667      : Dolan Fig.13 resolution
  - xmin = -750             : PML width 2.82*lambda (lambda=251 M), our >=2.5λ rule
  - eps_match = 0.002       : same stitch as the validated damped benchmark
  - T = 1e5                 : Dolan Fig.13 run length; Delta(ln|psi|) ~
                              Gamma*T/2 ~ 7.5e-3 over the fit window, >> the
                              ~1e-8 PML floor (floor test)
  - extraction = freqfilter : Dolan's own Fig.13 method; validated to 1.000x
                              Leaver on the a=0.9 damped case. logfit beats
                              (n=0/n=1 both grow), prony is Gamma-res-limited.

Success = Gamma comes out POSITIVE (growth) and near +1.5e-7. Full series saved
so logfit/prony/freqfilter can all be re-applied offline and a longer T=1e5
refinement can reuse the analysis.
"""
import numpy as np, sys, time
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
import evolve as ev

alpha, a = 0.42, 0.99
w = complex(leaver.solve(alpha, a))
print(f"Leaver target: E={w.real:.6f}  Gamma={w.imag:+.5e}  (want SIGN=+)", flush=True)

t0 = time.time()
o = ev.run(alpha=alpha, a=a, N=6000, xmin=-750., xmax=250., T=100000.,
           eps_match=0.002, gmax=1.5, nrec=16000)
print(f"run done in {time.time()-t0:.0f}s  probe r={o['probe_r']:.2f}  dx={o['dx']:.4f}"
      f"  x_match={o['x_match']:.1f}  PML={o['x_match']-(-750):.0f}={(o['x_match']-(-750))/251.3:.2f}lambda", flush=True)

t, psi = o["t"], o["psi"]
np.savez("stage_d_full.npz", t=t, psi=psi, w_leaver=w,
         phi_flux=o["phi_flux"], phidot_flux=o["phidot_flux"],
         phix_flux=o["phix_flux"], flux_r=o["flux_r"])
print("saved stage_d_full.npz", flush=True)

print(f"\n{'method':<26}{'E':>12}{'Gamma':>14}{'ratio':>10}", flush=True)
print(f"{'Leaver':<26}{w.real:>12.6f}{w.imag:>14.4e}{'1.000x':>10}", flush=True)
E1, G1 = ev.logfit(t, psi)
print(f"{'logfit':<26}{E1:>12.6f}{G1:>14.4e}{G1/w.imag:>9.3f}x", flush=True)
try:
    E2, G2 = ev.prony(t, psi, anchor=w)
    print(f"{'prony (anchored)':<26}{E2:>12.6f}{G2:>14.4e}{G2/w.imag:>9.3f}x", flush=True)
except Exception as e:
    print(f"prony failed: {e}", flush=True)
for kb in [3.0, 5.0, 8.0]:
    E3, G3 = ev.freqfilter(t, psi, E_peak=w.real, kbins=kb)
    print(f"{'freqfilter kbins='+str(kb):<26}{E3:>12.6f}{G3:>14.4e}{G3/w.imag:>9.3f}x", flush=True)
