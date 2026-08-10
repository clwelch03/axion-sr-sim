#!/usr/bin/env python
"""Exact-answer sanity tests for the coupled (l, m=1) tower (evolve_tower.py).

Four falsifiable predictions, none depending on fit quality:

  [1] REGRESSION: with l_max=1 (one channel) the tower must reproduce the
      single-mode evolve.run to MACHINE PRECISION -- same psi(t), same
      cloud_norm(t) -- for both lambda=0 and lambda!=0. This validates that the
      refactor changed nothing about the physics of the existing mode
      (Qbar^1_{111}=1 and the operator is the same single-mode operator).

  [2] SOURCING: only l=1 is seeded, and the off-diagonal LINEAR coupling is
      dropped, so l=3 can be populated ONLY by the cubic. Hence:
        (a) at lambda=0, ||psi_3|| is EXACTLY 0 for all time;
        (b) at lambda!=0, ||psi_3|| grows, LINEAR in lambda (a cubic source
            ~ lambda*psi_1^3): ||psi_3||(2L)/||psi_3||(L) -> 2 while l=3 is passive.

  [3] TOWER HERMITICITY: the multi-channel cubic is norm-conserving (angular.py
      shows sum_j Im(psi_j* F_j)=0), so it REDISTRIBUTES norm among channels but
      does no net work -- the TOTAL interior-norm growth rate is unchanged by
      lambda. Gamma_total(lambda) == Gamma_total(0) (the tower analog of the
      single-mode Hermiticity test).

  [4] CASCADE selection rule: with l_max=5 and only l=1 seeded, l=5 is
      SECOND-GENERATION -- three l=1's reach only l<=3 (Qbar^5_{111}=0), so l=5
      is sourced only after l=3 fills. Its amplitude relative to l=3 must GROW in
      time (l=5 catching up), not be a fixed ratio.

Cheap decisive tests (damped a=0.9, alpha=0.5, short T), not production runs.
"""
import numpy as np, argparse, sys, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))      # radial1d/
sys.path.insert(0, _HERE)                        # nonlinear/
import evolve
import evolve_tower


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.50)
    ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--T", type=float, default=200.)
    ap.add_argument("--N", type=int, default=1600)
    ap.add_argument("--xmin", type=float, default=-45.)
    ap.add_argument("--xmax", type=float, default=100.)
    ap.add_argument("--gmax", type=float, default=1.0)
    ap.add_argument("--eps", type=float, default=0.02)
    ap.add_argument("--lambda_ref", type=float, default=0.02)
    ap.add_argument("--amp", type=float, default=3.0)
    ap.add_argument("--plot", type=str, default=None)
    q = ap.parse_args()
    grid = dict(xmin=q.xmin, xmax=q.xmax, N=q.N, T=q.T, gmax=q.gmax, eps_match=q.eps)
    L = q.lambda_ref
    print(f"# tower sanity test  a={q.a} alpha={q.alpha} T={q.T} N={q.N} amp={q.amp} lambda_ref={L}")
    results = []

    # ---- [1] REGRESSION: l_max=1 tower == single-mode evolve.run (machine precision) ----
    print("\n[1] REGRESSION (l_max=1 tower == evolve.run):")
    reg_ok = True
    for lam in (0.0, L):
        single = evolve.run(q.alpha, q.a, q.xmin, q.xmax, q.N, q.T, gmax=q.gmax,
                            eps_match=q.eps, lambda_quartic=lam, init_amplitude=q.amp)
        tower = evolve_tower.run(q.alpha, q.a, 1, q.xmin, q.xmax, q.N, q.T, gmax=q.gmax,
                                eps_match=q.eps, lambda_quartic=lam, init_amplitude=q.amp)
        d_psi = np.max(np.abs(tower["psi"][0] - single["psi"]))
        d_cloud = np.max(np.abs(tower["cloud_norm"][0] - single["cloud_norm"]))
        scale = np.max(np.abs(single["psi"])) + 1e-300
        ok = d_psi < 1e-10 * scale and d_cloud < 1e-10 * (np.max(single["cloud_norm"]) + 1e-300)
        reg_ok &= ok
        print(f"    lambda={lam:+.3g}: max|dpsi|={d_psi:.2e}  max|dcloud|={d_cloud:.2e}  "
              f"{'PASS' if ok else 'FAIL'}")
    results.append(("[1] REGRESSION", reg_ok))

    # ---- [2] SOURCING: l=3 zero at lambda=0; linear in lambda otherwise ----
    print("\n[2] SOURCING (l=1 -> l=3 via the cubic):")
    t3_lin = evolve_tower.run(q.alpha, q.a, 3, **grid, lambda_quartic=0.0, init_amplitude=q.amp)
    l3_at_zero = np.max(np.abs(t3_lin["cloud_norm"][1]))
    zero_ok = l3_at_zero < 1e-14
    print(f"    lambda=0 : max||psi_3|| = {l3_at_zero:.2e}  (want ~0: nothing sources l=3)  "
          f"{'PASS' if zero_ok else 'FAIL'}")
    t3_L = evolve_tower.run(q.alpha, q.a, 3, **grid, lambda_quartic=L, init_amplitude=q.amp)
    t3_2L = evolve_tower.run(q.alpha, q.a, 3, **grid, lambda_quartic=2 * L, init_amplitude=q.amp)
    # ratio at an early record where l=3 is still passive (linear response to the l=1 source)
    k_early = len(t3_L["t"]) // 4
    n3_L = t3_L["cloud_norm"][1, k_early]; n3_2L = t3_2L["cloud_norm"][1, k_early]
    lam_ratio = n3_2L / n3_L if n3_L else np.nan
    lin_ok = abs(lam_ratio - 2.0) < 0.15
    print(f"    ||psi_3||(2L)/||psi_3||(L) @ t={t3_L['t'][k_early]:.0f} = {lam_ratio:.3f}  "
          f"(expect ~2, linear-in-lambda source)  {'PASS' if lin_ok else 'FAIL'}")
    results.append(("[2] SOURCING", zero_ok and lin_ok))

    # ---- [3] TOWER HERMITICITY: total-norm Gamma unchanged by lambda ----
    print("\n[3] TOWER HERMITICITY (total interior Gamma independent of lambda):")
    def total_gamma(run_output):
        total = np.sqrt(np.sum(run_output["cloud_norm"] ** 2, axis=0))   # sqrt(sum_l ||psi_l||^2)
        t_g, g_t = evolve.inst_gamma(run_output["t"], np.maximum(total, 1e-300))
        return np.mean(g_t[t_g > 0.6 * t_g[-1]])
    g_lin = total_gamma(t3_lin)          # lambda=0 (reuse from test 2)
    g_L = total_gamma(t3_L)              # lambda=L
    gamma_shift = abs(g_L - g_lin)
    herm_ok = gamma_shift < 0.03 * abs(g_lin)
    print(f"    Gamma_total(0)={g_lin:+.4e}  Gamma_total(L)={g_L:+.4e}  "
          f"|dGamma|={gamma_shift:.2e} (< 3% of |Gamma|={0.03*abs(g_lin):.2e}?)  "
          f"{'PASS' if herm_ok else 'FAIL'}")
    results.append(("[3] TOWER HERMITICITY", herm_ok))

    # ---- [4] CASCADE: l=5 is second-generation (grows relative to l=3 over time) ----
    print("\n[4] CASCADE selection rule (l=5 second-generation, l_max=5):")
    t5 = evolve_tower.run(q.alpha, q.a, 5, **grid, lambda_quartic=L, init_amplitude=q.amp)
    n3 = t5["cloud_norm"][1]; n5 = t5["cloud_norm"][2]
    k_mid = len(t5["t"]) // 2
    ratio_mid = n5[k_mid] / n3[k_mid] if n3[k_mid] else np.nan
    ratio_end = n5[-1] / n3[-1] if n3[-1] else np.nan
    cascade_ok = ratio_end > ratio_mid > 0
    print(f"    ||psi_5||/||psi_3|| : mid={ratio_mid:.3e}  end={ratio_end:.3e}  "
          f"(end>mid => l=5 fills after l=3?)  {'PASS' if cascade_ok else 'FAIL'}")
    results.append(("[4] CASCADE", cascade_ok))

    all_ok = all(ok for _, ok in results)
    print("\n  " + "  ".join(f"{name}:{'OK' if ok else 'X'}" for name, ok in results))
    print(f"  OVERALL: {'ALL PASS' if all_ok else 'CHECK FAILURES ABOVE'}")

    if q.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        ls = t5["ls"]
        for i, l in enumerate(ls):
            axes[0].plot(t5["t"], np.abs(t5["cloud_norm"][i]), label=f"l={l}")
        axes[0].set_yscale("log"); axes[0].set_xlabel("t / M")
        axes[0].set_ylabel("||psi_l|| (interior L2)"); axes[0].legend()
        axes[0].set_title(f"tower populations (lambda={L}, l_max=5)")
        axes[1].plot(t5["t"], n5 / np.maximum(n3, 1e-300))
        axes[1].set_xlabel("t / M"); axes[1].set_ylabel("||psi_5|| / ||psi_3||")
        axes[1].set_title("l=5 second-generation growth")
        fig.tight_layout(); fig.savefig(q.plot, dpi=110); print("saved", q.plot)


if __name__ == "__main__":
    main()
