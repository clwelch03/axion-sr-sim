#!/usr/bin/env python
"""One growth-peak tower run for the Gamma(A, l_max) multi-mode growth sweep.

Measures how the coupled (l, m=1) tower modifies the SUPERRADIANT GROWTH rate
of the l=1 mode as a function of seed amplitude A and truncation l_max, at the
Dolan/Leaver growth peak (a=0.99, alpha=0.42, Leaver Gamma=+1.5183e-7).

Method = the validated Stage-D pipeline (docs/TIER1_GAMMA_NLM_RESULTS.md 4.1):
Dolan near-horizon PML, T=1e5 M so Delta(ln|psi|)~Gamma*T/2 >> the ~1e-8 floor,
freqfilter extraction (Dolan Fig.13). Two changes for the tower:
  * seed = the CLEAN l=1 eigenmode (radial1d.eigenmode_seed) instead of a
    Gaussian -- suppresses the n=1-overtone beat that freqfilter was built to
    fight, so the per-channel rate is cleaner.
  * evolve_tower with the RWA quartic ON: the conservative cubic pumps the
    (growing) l=1 into the (strongly LINEARLY DAMPED) l=3, l=5 channels, which
    bleed down the horizon -- a NET growth-suppression that the single-mode
    (Hermitian) cubic cannot produce. l_max=1 is therefore the CONTROL: its
    growth rate must stay flat at +1.5e-7 (Hermitian => detuning only, no
    Gamma change); l_max>=3 is where Gamma(A) falls.

Headline number = the l=1 freqfilter Gamma. Differencing against the lambda=0
baseline (same grid) cancels the freqfilter/PML systematic (the ~1.056x
overshoot) to give a clean Delta_Gamma(A, l_max).
"""
import numpy as np, sys, os, time, argparse

_LR = "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates"
_RD = "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d"
sys.path.insert(0, _LR); sys.path.insert(0, _RD); sys.path.insert(0, os.path.join(_RD, "nonlinear"))
import leaver
import evolve as ev
import evolve_tower as et

# growth peak (Stage-D definitive config)
ALPHA, A_SPIN = 0.42, 0.99
XMIN, XMAX, EPS_MATCH, GMAX = -750., 250., 0.002, 1.5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l_max", type=int, required=True, help="tower truncation (1=control, 3, 5, ...)")
    ap.add_argument("--amp", type=float, required=True, help="l=1 eigenmode seed amplitude A")
    ap.add_argument("--lam", type=float, required=True, help="lambda_quartic (0 = linear baseline)")
    ap.add_argument("--T", type=float, default=1e5, help="evolution time in M (Dolan Fig.13 = 1e5)")
    ap.add_argument("--N", type=int, default=6000)
    ap.add_argument("--nrec", type=int, default=16000)
    ap.add_argument("--out", required=True, help="output .npz path")
    q = ap.parse_args()

    w = complex(leaver.solve(ALPHA, A_SPIN))
    print(f"Leaver (l=1): E={w.real:.6f}  Gamma={w.imag:+.5e}   "
          f"[l_max={q.l_max} A={q.amp} lambda={q.lam} T={q.T:.0e} N={q.N}]", flush=True)

    t0 = time.time()
    o = et.run(ALPHA, A_SPIN, q.l_max, xmin=XMIN, xmax=XMAX, N=q.N, T=q.T,
               eps_match=EPS_MATCH, gmax=GMAX, nrec=q.nrec,
               lambda_quartic=q.lam, init_amplitude=q.amp, seed="eigenmode")
    dt_run = time.time() - t0
    pml_lambda = (o["x_match"] - XMIN) / 251.3
    print(f"run done in {dt_run:.0f}s  dx={o['dx']:.4f}  x_match={o['x_match']:.1f}  "
          f"PML={o['x_match']-XMIN:.0f} = {pml_lambda:.2f}lambda  "
          f"finite={np.all(np.isfinite(o['psi']))}", flush=True)

    t, psi, ls = o["t"], o["psi"], o["ls"]
    np.savez(q.out, t=t, psi=psi, cloud_norm=o["cloud_norm"], ls=ls, w_leaver=w,
             alpha=ALPHA, a=A_SPIN, l_max=q.l_max, amp=q.amp, lam=q.lam, T=q.T, N=q.N,
             dx=o["dx"], x_match=o["x_match"], run_seconds=dt_run)
    print(f"saved {q.out}", flush=True)

    # per-channel rate extraction (resonant tower => all channels share the l=1
    # frequency, so freqfilter is centred on the l=1 Leaver E for every channel).
    print(f"\n{'chan':<6}{'method':<20}{'E':>12}{'Gamma':>15}{'ratio(l1lin)':>14}", flush=True)
    print(f"{'--':<6}{'Leaver l=1':<20}{w.real:>12.6f}{w.imag:>15.4e}{'1.000x':>14}", flush=True)
    for i, l in enumerate(ls):
        p = psi[i]
        for kb in (3.0, 5.0, 8.0):
            E, G = ev.freqfilter(t, p, E_peak=w.real, kbins=kb)
            print(f"{('l='+str(l)):<6}{('freqfilter kb='+str(kb)):<20}"
                  f"{E:>12.6f}{G:>15.4e}{G/w.imag:>13.3f}x", flush=True)
        E1, G1 = ev.logfit(t, p)
        print(f"{('l='+str(l)):<6}{'logfit':<20}{E1:>12.6f}{G1:>15.4e}{G1/w.imag:>13.3f}x", flush=True)


if __name__ == "__main__":
    main()
