#!/usr/bin/env python
"""Exact-answer sanity tests for the RWA quartic self-interaction
(evolve.run lambda_quartic).

The added term -Delta*lambda_quartic*|psi|^2*psi is Hermitian, so it makes two
falsifiable predictions that DON'T depend on any fit quality:

  (1) HERMITICITY: it is norm-conserving, so it CANNOT change Gamma directly.
      At fixed amplitude, Gamma(lambda) must equal Gamma(0) (growth comes only
      from the horizon). Any Gamma shift = a non-Hermitian discretization bug
      (the check the README's working-style section asks for).

  (2) DETUNING SCALING: it shifts the real frequency by
      delta_E = E(lambda) - E(0), with delta_E linear in lambda (sign flips) and
      quadratic in amplitude (delta_E ~ lambda*<|psi|^2> ~ lambda*amp^2).

Run on the damped benchmark (a=0.9, alpha=0.5), short T so |psi| is ~stationary
over the fit window. This is a cheap decisive test, not a production run.
"""
import numpy as np, argparse, sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d")
from evolve import run, logfit, inst_gamma

def measure(alpha, a, lambda_quartic, init_amplitude, T, N, xmin, xmax, gmax, eps):
    """Evolve once and return (E, Gamma, run_output) via a single-mode log-fit."""
    run_output = run(alpha, a, xmin, xmax, N, T, gmax=gmax, eps_match=eps,
                     lambda_quartic=lambda_quartic, init_amplitude=init_amplitude)
    freq, gamma = logfit(run_output["t"], run_output["psi"], tmin_frac=0.35)
    return freq, gamma, run_output

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.50)
    ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--T", type=float, default=400.)
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--xmin", type=float, default=-70.)
    ap.add_argument("--xmax", type=float, default=220.)
    ap.add_argument("--gmax", type=float, default=1.0)
    ap.add_argument("--eps", type=float, default=0.02)
    ap.add_argument("--lambda_ref", type=float, default=0.02,
                    help="reference |lambda| for the scaling tests")
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--plot", type=str, default=None)
    q = ap.parse_args()
    grid_kw = dict(T=q.T, N=q.N, xmin=q.xmin, xmax=q.xmax, gmax=q.gmax, eps=q.eps)
    lambda_ref = q.lambda_ref

    print(f"# self-interaction sanity test  a={q.a} alpha={q.alpha} T={q.T} N={q.N}")
    # baseline (linear)
    freq_linear, gamma_linear, run_linear = measure(q.alpha, q.a, 0.0, q.amp, **grid_kw)
    omega_leaver = run_linear["w_ref"]
    print(f"  Leaver           : E={omega_leaver.real:.6f}  Gamma={omega_leaver.imag:+.4e}")
    print(f"  lambda=0 (base)  : E={freq_linear:.6f}  Gamma={gamma_linear:+.4e}")

    # --- Test 1: Hermiticity (Gamma unchanged) + Test 2a: dE sign flip ---
    freq_pos, gamma_pos, _ = measure(q.alpha, q.a, +lambda_ref, q.amp, **grid_kw)
    freq_neg, gamma_neg, _ = measure(q.alpha, q.a, -lambda_ref, q.amp, **grid_kw)
    print(f"  lambda=+{lambda_ref:<7g}  : E={freq_pos:.6f}  Gamma={gamma_pos:+.4e}   "
          f"dE={freq_pos-freq_linear:+.3e}")
    print(f"  lambda=-{lambda_ref:<7g}  : E={freq_neg:.6f}  Gamma={gamma_neg:+.4e}   "
          f"dE={freq_neg-freq_linear:+.3e}")

    gamma_shift = max(abs(gamma_pos-gamma_linear), abs(gamma_neg-gamma_linear))
    hermiticity_ok = gamma_shift < 0.03*abs(gamma_linear)   # Gamma stable to <3% of |Gamma|
    print(f"\n  [1] HERMITICITY : max|dGamma|={gamma_shift:.2e}  "
          f"(< 3% of |Gamma|={0.03*abs(gamma_linear):.2e}?)  "
          f"{'PASS' if hermiticity_ok else 'FAIL'}")

    dfreq_pos, dfreq_neg = freq_pos-freq_linear, freq_neg-freq_linear
    sign_flip_ok = (dfreq_pos*dfreq_neg < 0 and
                    abs(abs(dfreq_pos)-abs(dfreq_neg)) < 0.25*max(abs(dfreq_pos), abs(dfreq_neg)))
    print(f"  [2a] dE SIGN-FLIP: dE(+lambda)={dfreq_pos:+.3e}  dE(-lambda)={dfreq_neg:+.3e}  "
          f"(opposite & ~equal?)  {'PASS' if sign_flip_ok else 'FAIL'}")

    # --- Test 2b: dE ~ amp^2 at fixed lambda ---
    freq_amp1, _, _ = measure(q.alpha, q.a, +lambda_ref, q.amp,     **grid_kw)
    freq_amp2, _, _ = measure(q.alpha, q.a, +lambda_ref, 2.0*q.amp, **grid_kw)
    amp_ratio = (freq_amp2-freq_linear)/(freq_amp1-freq_linear) if (freq_amp1-freq_linear) else np.nan
    amp_scaling_ok = abs(amp_ratio-4.0) < 1.0        # expect ~4 (quadratic)
    print(f"  [2b] dE ~ amp^2  : dE(2amp)/dE(amp)={amp_ratio:.2f}  (expect ~4)  "
          f"{'PASS' if amp_scaling_ok else 'FAIL'}")

    # --- Test 2c: dE ~ lambda at fixed amp ---
    freq_half, _, _ = measure(q.alpha, q.a, +0.5*lambda_ref, q.amp, **grid_kw)
    lambda_ratio = (freq_pos-freq_linear)/(freq_half-freq_linear) if (freq_half-freq_linear) else np.nan
    lambda_scaling_ok = abs(lambda_ratio-2.0) < 0.5  # expect ~2 (linear)
    print(f"  [2c] dE ~ lambda : dE(lambda)/dE(lambda/2)={lambda_ratio:.2f}  (expect ~2)  "
          f"{'PASS' if lambda_scaling_ok else 'FAIL'}")

    all_ok = hermiticity_ok and sign_flip_ok and amp_scaling_ok and lambda_scaling_ok
    print(f"\n  OVERALL: {'ALL PASS' if all_ok else 'CHECK FAILURES ABOVE'}")

    if q.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        t_gamma_lin, gamma_t_lin = inst_gamma(run_linear["t"], run_linear["cloud_norm"])
        run_pos = run(q.alpha, q.a, q.xmin, q.xmax, q.N, q.T, gmax=q.gmax,
                      eps_match=q.eps, lambda_quartic=+lambda_ref, init_amplitude=q.amp)
        t_gamma_pos, gamma_t_pos = inst_gamma(run_pos["t"], run_pos["cloud_norm"])
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].axhline(omega_leaver.imag, ls="--", c="k", lw=1, label="Leaver Gamma")
        axes[0].plot(t_gamma_lin, gamma_t_lin, label="lambda=0")
        axes[0].plot(t_gamma_pos, gamma_t_pos, label=f"lambda=+{lambda_ref}")
        axes[0].set_xlabel("t / M"); axes[0].set_ylabel("Gamma(t) = d ln||psi||/dt")
        axes[0].legend()
        axes[1].plot(run_linear["t"], np.abs(run_linear["cloud_norm"]), label="lambda=0")
        axes[1].plot(run_pos["t"], np.abs(run_pos["cloud_norm"]), label=f"lambda=+{lambda_ref}")
        axes[1].set_yscale("log"); axes[1].set_xlabel("t / M"); axes[1].set_ylabel("||psi||")
        axes[1].legend()
        fig.tight_layout(); fig.savefig(q.plot, dpi=110); print("saved", q.plot)

if __name__ == "__main__":
    main()
