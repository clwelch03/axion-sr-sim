#!/usr/bin/env python
"""Detuning sweep: calibrate the RWA-quartic frequency shift and PREDICT the
detuning-driven saturation amplitude -- the cheap-decisive deliverable of the
single-mode nonlinear step (README Sec 10.5 [A], item 2).

Physics: the Hermitian cubic shifts the real frequency by
    delta_E(A) = k_detune * lambda * A^2        (perturbative RWA)
and does NOT change Gamma except through that shift. Superradiance runs while
E < m*OmH and shuts off at E = m*OmH, so a growing mode self-saturates (via
detuning ALONE) at the amplitude where

    delta_E(A_sat) = m*OmH - E0        (the gap the mode must close).

k_detune is a property of the mode's radial profile: we MEASURE it cheaply on a
short run (differencing lambda=+/-L against lambda=0 at fixed amplitude cancels
the absolute-E fit systematic), verify the A^2 * lambda law over a decade of
amplitude, then report the predicted saturation gap. This lets us decide whether
a full run-to-saturation is even feasible before paying for it.

CAVEAT: for a GROWING mode E0 < m*OmH; lambda>0 (repulsive) RAISES E toward the
threshold => self-limiting. lambda<0 (attractive/axion) LOWERS E, away from
threshold => detuning does NOT saturate it (its real saturation is collapse /
multi-mode, not captured here -- README Sec 10.1).
"""
import numpy as np, argparse, sys
# evolve.py lives in the parent radial1d/ dir (this script sits in nonlinear/);
# importing evolve transitively pulls in leaver via evolve's own path insert.
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
from evolve import run, logfit, horizons

def frequency_at(alpha, a, lambda_quartic, init_amplitude, grid_kw):
    """Real frequency E from a single short evolution (log-fit phase slope)."""
    run_output = run(alpha, a, lambda_quartic=lambda_quartic,
                     init_amplitude=init_amplitude, **grid_kw)
    freq, _ = logfit(run_output["t"], run_output["psi"], tmin_frac=0.35)
    return freq, run_output

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.50)
    ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--T", type=float, default=400.)
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--xmin", type=float, default=-70.)
    ap.add_argument("--xmax", type=float, default=220.)
    ap.add_argument("--gmax", type=float, default=1.0)
    ap.add_argument("--eps_match", type=float, default=0.02)
    ap.add_argument("--lambda_ref", type=float, default=0.02,
                    help="|lambda| used to probe the shift (differenced +/- )")
    ap.add_argument("--amplitudes", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 4.0])
    ap.add_argument("--plot", type=str, default=None)
    q = ap.parse_args()
    grid_kw = dict(T=q.T, N=q.N, xmin=q.xmin, xmax=q.xmax, gmax=q.gmax,
                   eps_match=q.eps_match)
    lambda_ref = q.lambda_ref

    r_plus, _ = horizons(q.a); horizon_freq = 1 * q.a/(2*r_plus)   # m=1
    _, run0 = frequency_at(q.alpha, q.a, 0.0, q.amplitudes[0], grid_kw)
    omega_leaver = run0["w_ref"]
    E0 = omega_leaver.real
    gap = horizon_freq - E0            # >0 for a growing (superradiant) mode
    print(f"# detuning sweep  a={q.a} alpha={q.alpha} T={q.T} N={q.N} |lambda|={lambda_ref}")
    print(f"  Leaver E0={E0:.6f}  Gamma={omega_leaver.imag:+.3e}  "
          f"m*OmH={horizon_freq:.6f}  gap(m*OmH-E0)={gap:+.4e}  "
          f"({'GROWING' if gap>0 else 'DAMPED'})")
    print(f"  {'amp':>6} {'amp^2':>8} {'dE(+lam)':>12} {'dE(-lam)':>12} "
          f"{'asym':>9} {'k_detune':>11}")

    amp_sq, dE_pos_list, dE_neg_list = [], [], []
    for amplitude in q.amplitudes:
        E_base, _ = frequency_at(q.alpha, q.a, 0.0,        amplitude, grid_kw)
        E_pos,  _ = frequency_at(q.alpha, q.a, +lambda_ref, amplitude, grid_kw)
        E_neg,  _ = frequency_at(q.alpha, q.a, -lambda_ref, amplitude, grid_kw)
        dE_pos, dE_neg = E_pos - E_base, E_neg - E_base
        # k from the antisymmetric part (cancels any residual even-in-lambda drift)
        k_detune = 0.5*(dE_pos - dE_neg)/(lambda_ref*amplitude**2)
        asymmetry = abs(dE_pos + dE_neg)/max(abs(dE_pos), abs(dE_neg), 1e-300)
        amp_sq.append(amplitude**2); dE_pos_list.append(dE_pos); dE_neg_list.append(dE_neg)
        print(f"  {amplitude:6.3f} {amplitude**2:8.3f} {dE_pos:+12.3e} "
              f"{dE_neg:+12.3e} {asymmetry:9.1e} {k_detune:+11.3e}")

    amp_sq = np.array(amp_sq)
    # slope of the antisymmetric shift vs amp^2 -> single calibrated k_detune
    dE_anti = 0.5*(np.array(dE_pos_list) - np.array(dE_neg_list))
    slope, intercept = np.polyfit(amp_sq, dE_anti, 1)
    k_detune_fit = slope/lambda_ref
    residual = dE_anti - (slope*amp_sq + intercept)
    linearity = 1 - np.sum(residual**2)/np.sum((dE_anti-dE_anti.mean())**2)
    print(f"\n  FIT dE_anti = slope*amp^2 + b : slope={slope:+.4e}  "
          f"k_detune={k_detune_fit:+.4e}  R^2(amp^2 law)={linearity:.5f}")

    # saturation-amplitude prediction: ONLY meaningful for a GROWING mode
    # (gap>0). A damped mode (gap<0) decays -- it never grows to saturation, so
    # detuning saturation does not apply.
    if gap <= 0:
        print(f"\n  PREDICTION: N/A -- this mode is DAMPED (gap<=0), it does not grow "
              f"to a saturation amplitude. Run at a growing config (E0<m*OmH, e.g. "
              f"a=0.99 alpha~0.42) for a physical A_sat.")
    elif k_detune_fit != 0:
        # detuning closes the gap when k*lambda*A^2 = gap; needs sign(lambda)=sign(gap/k)
        lam_sign = np.sign(gap/k_detune_fit)
        A_sat_sq = gap/(k_detune_fit*lam_sign*abs(lambda_ref))
        print(f"\n  PREDICTION (detuning-only saturation): closing delta_E={gap:+.4e} "
              f"needs lambda={lam_sign:+.0f}*|lambda| (repulsive raises E toward m*OmH); "
              f"A_sat^2 = {A_sat_sq:.3e} at |lambda|={abs(lambda_ref)}, "
              f"i.e. A_sat = {np.sqrt(A_sat_sq):.3e}*sqrt({abs(lambda_ref)}/|lambda|).")
        print(f"  (The OTHER sign detunes AWAY from threshold and does NOT saturate "
              f"by detuning -- its saturation is collapse/multi-mode, not modeled here.)")

    if q.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.plot(amp_sq, np.array(dE_pos_list), "o-", label=f"lambda=+{lambda_ref}")
        ax.plot(amp_sq, np.array(dE_neg_list), "s-", label=f"lambda=-{lambda_ref}")
        ax.plot(amp_sq, slope*amp_sq+intercept, "k--", lw=1,
                label=f"fit k={k_detune_fit:.2e}")
        ax.set_xlabel("init_amplitude^2"); ax.set_ylabel("delta_E = E(lambda)-E(0)")
        ax.set_title(f"RWA detuning  a={q.a} alpha={q.alpha}"); ax.legend()
        ax.axhline(0, c="0.7", lw=0.8)
        fig.tight_layout(); fig.savefig(q.plot, dpi=110); print("saved", q.plot)

if __name__ == "__main__":
    main()
