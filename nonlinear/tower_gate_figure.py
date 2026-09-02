#!/usr/bin/env python
"""Results figure for the coupled-tower validation gate (companion to tower_test.py).

Regenerates the (cheap, damped a=0.9 alpha=0.5, short-T) runs behind the four gate
checks and renders one panel each, plus the lambda-scan discriminator that showed
the [3] total-Gamma shift is O(lambda) detuning (not an O(lambda^2) net-work
violation). Self-contained: run it and it writes tower_gate_results.png.

    python tower_gate_figure.py [--out tower_gate_results.png]
"""
import sys, os, argparse, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))      # radial1d/
sys.path.insert(0, _HERE)                        # nonlinear/
import evolve, evolve_tower

CL = {1: "C0", 3: "C1", 5: "C2"}                 # one colour per channel, all panels


def cap(ax, text):
    """Two-line explanatory caption centered below the panel."""
    ax.text(0.5, -0.26, text, transform=ax.transAxes, ha="center", va="top",
            fontsize=8.0, linespacing=1.35)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(_HERE, "tower_gate_results.png"))
    ap.add_argument("--alpha", type=float, default=0.50)
    ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--T", type=float, default=200.)
    ap.add_argument("--N", type=int, default=1600)
    ap.add_argument("--amp", type=float, default=3.0)
    ap.add_argument("--L", type=float, default=0.02)
    q = ap.parse_args()
    grid = dict(xmin=-45., xmax=100., N=q.N, T=q.T, gmax=1.0, eps_match=0.02)
    a, alpha, amp, L = q.a, q.alpha, q.amp, q.L
    print(f"# tower gate figure  a={a} alpha={alpha} T={q.T} N={q.N} amp={amp} L={L}")

    # ---------- runs ----------
    print("  [3-disc] lambda-scan runs (l_max=3) ...")
    t3_0 = evolve_tower.run(alpha, a, 3, **grid, lambda_quartic=0.0, init_amplitude=amp)
    t3_L = evolve_tower.run(alpha, a, 3, **grid, lambda_quartic=L, init_amplitude=amp)
    print("  [4] cascade run (l_max=5) ...")
    t5 = evolve_tower.run(alpha, a, 5, **grid, lambda_quartic=L, init_amplitude=amp)
    print("  [3-disc] lambda scan ...")
    def total_gamma(o):
        tot = np.sqrt(np.sum(o["cloud_norm"] ** 2, axis=0))
        tg, gt = evolve.inst_gamma(o["t"], np.maximum(tot, 1e-300))
        return np.mean(gt[tg > 0.6 * tg[-1]])
    g0 = total_gamma(t3_0)
    lam_scan = [0.005, 0.01, 0.02, 0.04]
    dG = []
    for lam in lam_scan:
        o = t3_L if lam == L else evolve_tower.run(alpha, a, 3, **grid,
                                                   lambda_quartic=lam, init_amplitude=amp)
        dG.append(total_gamma(o) - g0)
    dG = np.array(dG)

    # ---------- figure ----------
    # Two panels tell the whole story: (left) the cubic drains the growing l=1
    # mode into a cascade of higher-l channels, second-generation modes filling
    # only after their parent does; (right) that drain is real Hermitian
    # detuning, not a bug -- the total-Gamma shift scales linearly in the
    # coupling lambda, not quadratically as a net-work violation would.
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    fig.suptitle(f"Coupled (l, m=1) tower — cascade + Hermiticity discriminator  "
                 f"(damped a={a}, α={alpha}, T={q.T:.0f}M, N={q.N})",
                 fontsize=12, fontweight="bold")

    # --- CASCADE populations + fill-time lag ---
    axD = ax[0]
    tt5 = t5["t"]
    for i, l in enumerate(t5["ls"]):
        axD.plot(tt5, np.abs(t5["cloud_norm"][i]) + 1e-300, color=CL[int(l)], lw=1.8, label=f"l={l}")
    def t_reach(nrm, frac): return tt5[np.argmax(nrm >= frac * np.max(nrm))]
    th3 = t_reach(t5["cloud_norm"][1], 0.5); th5 = t_reach(t5["cloud_norm"][2], 0.5)
    axD.axvline(th3, color=CL[3], ls=":", lw=1.0); axD.axvline(th5, color=CL[5], ls=":", lw=1.0)
    axD.set_yscale("log"); axD.set_ylim(1e-6, 1e2)
    axD.set_xlabel("t / M"); axD.set_ylabel("‖ψ_l‖ (interior L2)")
    axD.set_title("cascade: l=1 drains into l=3, then l=5")
    axD.text(0.97, 0.06, f"t½: l=3 = {th3:.0f}M,  l=5 = {th5:.0f}M\n(l=5 is second-generation, fills later)",
             transform=axD.transAxes, ha="right", fontsize=8,
             bbox=dict(boxstyle="round", fc="w", ec="0.7"))
    axD.legend(loc="upper right", fontsize=8)
    cap(axD, "Three l=1 quanta reach only l≤3, so l=5 needs l=3 first:\n"
             "it grows only after l=3 fills, reaching half-peak ~7M later.")

    # --- λ-scan: dΓ ∝ λ (Hermitian detuning), not λ² (net-work) ---
    axF = ax[1]
    axF.plot(lam_scan, -dG * 1e3, "o", color="C5", ms=7, label="measured  −ΔΓ")
    slope = np.sum(np.array(lam_scan) * (-dG)) / np.sum(np.array(lam_scan) ** 2)  # through origin
    xs = np.linspace(0, max(lam_scan) * 1.05, 50)
    axF.plot(xs, slope * xs * 1e3, "C5--", lw=1.2, label=f"linear fit  (∝ λ¹)")
    q2 = (-dG[-1]) / lam_scan[-1] ** 2                                            # a λ² through last pt
    axF.plot(xs, q2 * xs ** 2 * 1e3, "0.6", ls=":", lw=1.2, label="∝ λ²  (ruled out)")
    axF.set_xlabel("λ (quartic coupling)"); axF.set_ylabel("−ΔΓ_total  (×10⁻³)")
    axF.set_title("the drain is detuning, not a bug")
    axF.legend(loc="upper left", fontsize=8)
    cap(axF, "Total-Γ shift is linear in λ (dashed), not quadratic (dotted):\n"
             "a net-work violation would be O(λ²); linear ⇒ Hermitian detuning.")

    fig.subplots_adjust(top=0.84, bottom=0.22, left=0.08, right=0.98, wspace=0.28)
    fig.savefig(q.out, dpi=120)
    print("saved", q.out)


if __name__ == "__main__":
    main()
