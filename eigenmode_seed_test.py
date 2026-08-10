#!/usr/bin/env python
"""Validate the clean eigenmode seed (radial1d.eigenmode_seed) against the Gaussian.

Exact-answer test (CLAUDE.md working style): on a LINEAR (lambda=0) run a correct
eigenmode seed must reproduce the Leaver growth rate from early times with no
shedding transient -- the Gaussian seed rings as it settles onto the true mode.
We run both seeds at lambda=0 on the damped config (a=0.9, alpha=0.5, Leaver
Gamma=-5.2e-4) and compare (i) the late-time rate vs Leaver and (ii) the transient
wobble std(inst_gamma)/|Gamma|. A single log-amplitude fit of the (smooth)
eigenmode cloud_norm recovers Gamma to high accuracy.

    python eigenmode_seed_test.py [--plot eigenmode_seed_test.png]
"""
import sys, os, argparse, numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import evolve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.50)
    ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--T", type=float, default=400.)
    ap.add_argument("--N", type=int, default=1600)
    ap.add_argument("--plot", default=None)
    q = ap.parse_args()
    grid = dict(xmin=-45., xmax=100., N=q.N, T=q.T, gmax=1.0, eps_match=0.02)
    print(f"# eigenmode-seed test  a={q.a} alpha={q.alpha} T={q.T} N={q.N}  (lambda=0)")
    gauss = evolve.run(q.alpha, q.a, **grid, lambda_quartic=0.0, init_amplitude=1.0, seed="gaussian")
    eig = evolve.run(q.alpha, q.a, **grid, lambda_quartic=0.0, init_amplitude=1.0, seed="eigenmode")
    GL = gauss["w_ref"].imag
    print(f"  Leaver Gamma = {GL:+.4e}")
    results = {}
    for name, o in (("gaussian", gauss), ("eigenmode", eig)):
        tg, gt = evolve.inst_gamma(o["t"], np.maximum(o["cloud_norm"], 1e-300))
        late = tg > 0.6 * tg[-1]
        g_late = np.mean(gt[late])
        wobble = np.std(gt[tg > 10]) / abs(GL)
        # single log-amplitude fit over the settled second half (phase-free rate)
        mfit = o["t"] > 0.5 * o["t"][-1]
        g_fit = np.polyfit(o["t"][mfit], np.log(o["cloud_norm"][mfit]), 1)[0]
        results[name] = (g_late, g_fit, wobble)
        print(f"  {name:9s}: Gamma_late={g_late:+.4e} ({g_late/GL:+.2f}x)  "
              f"Gamma_logfit={g_fit:+.4e} ({g_fit/GL:+.2f}x)  wobble={wobble:.2f}")
    ok = (abs(results["eigenmode"][1] / GL - 1) < 0.1
          and results["eigenmode"][2] < 0.1 * results["gaussian"][2])
    print(f"  -> eigenmode logfit within 10% of Leaver AND >10x less wobble: "
          f"{'PASS' if ok else 'CHECK'}")

    if q.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for name, o, c in (("gaussian", gauss, "C0"), ("eigenmode", eig, "C1")):
            ax[0].plot(o["t"], o["cloud_norm"], c, label=name)
            tg, gt = evolve.inst_gamma(o["t"], np.maximum(o["cloud_norm"], 1e-300))
            ax[1].plot(tg, gt, c, label=name)
        ax[0].set_yscale("log"); ax[0].set_xlabel("t / M"); ax[0].set_ylabel("cloud_norm")
        ax[0].legend(); ax[0].set_title(f"interior L2 norm (λ=0, a={q.a}, α={q.alpha})")
        ax[1].axhline(GL, color="k", ls="--", lw=1, label=f"Leaver {GL:.2e}")
        ax[1].set_xlabel("t / M"); ax[1].set_ylabel("inst_gamma")
        ax[1].set_ylim(GL * 3, abs(GL) * 3); ax[1].legend()
        ax[1].set_title("instantaneous growth rate")
        fig.tight_layout(); fig.savefig(q.plot, dpi=110); print("saved", q.plot)


if __name__ == "__main__":
    main()
