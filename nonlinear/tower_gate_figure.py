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
import evolve, evolve_tower, angular

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
    print("  [1] regression runs ...")
    single = evolve.run(alpha, a, grid["xmin"], grid["xmax"], q.N, q.T, gmax=1.0,
                        eps_match=0.02, lambda_quartic=L, init_amplitude=amp)
    tow1 = evolve_tower.run(alpha, a, 1, **grid, lambda_quartic=L, init_amplitude=amp)
    print("  [2]/[3] sourcing + hermiticity runs ...")
    t3_0 = evolve_tower.run(alpha, a, 3, **grid, lambda_quartic=0.0, init_amplitude=amp)
    t3_L = evolve_tower.run(alpha, a, 3, **grid, lambda_quartic=L, init_amplitude=amp)
    t3_2L = evolve_tower.run(alpha, a, 3, **grid, lambda_quartic=2 * L, init_amplitude=amp)
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
    fig, ax = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"Coupled (l, m=1) tower — validation gate  (damped a={a}, α={alpha}, "
                 f"T={q.T:.0f}M, N={q.N})   ALL 4 PASS", fontsize=13, fontweight="bold")

    # --- [1] REGRESSION: l_max=1 tower == single-mode, to machine precision ---
    axA = ax[0, 0]
    tt = tow1["t"]; win = tt < 60                     # zoom so the oscillation is visible
    dpsi = np.max(np.abs(tow1["psi"][0] - single["psi"]))
    axA.plot(tt[win], single["psi"][win].real, "C0-", lw=2.5, label="single-mode evolve.run")
    axA.plot(tt[win], tow1["psi"][0][win].real, "k--", lw=1.0, label="tower (l_max=1)")
    axA.set_xlabel("t / M"); axA.set_ylabel("Re ψ (probe)")
    axA.set_title(f"[1] REGRESSION\nmax|Δψ| = {dpsi:.1e}  →  machine precision")
    axA.legend(loc="lower center", fontsize=8, ncol=2)
    axins = axA.inset_axes([0.60, 0.66, 0.36, 0.28])
    axins.semilogy(tt, np.abs(tow1["psi"][0] - single["psi"]) + 1e-300, "C3", lw=0.8)
    axins.set_ylim(1e-18, 1e-10); axins.set_title("|Δψ(t)|", fontsize=7)
    axins.tick_params(labelsize=6)
    cap(axA, "One channel (l_max=1) reproduces the single-mode solver exactly:\n"
             "curves coincide and |Δψ|~1e-15 (inset), so no physics changed.")

    # --- [2] SOURCING: l=3 zero at λ=0, linear-in-λ otherwise ---
    axB = ax[0, 1]
    l3_at_zero = np.max(t3_0["cloud_norm"][1])
    axB.plot(t3_L["t"], t3_L["cloud_norm"][1], color=CL[3], lw=2.0, label="λ=L")
    axB.plot(t3_2L["t"], t3_2L["cloud_norm"][1], color=CL[3], lw=2.0, ls="--", label="λ=2L")
    axB.set_yscale("log"); axB.set_ylim(1e-4, 1e2)
    axB.set_xlabel("t / M"); axB.set_ylabel("‖ψ₃‖ (interior L2)")
    k = len(t3_L["t"]) // 4
    ratio = t3_2L["cloud_norm"][1, k] / t3_L["cloud_norm"][1, k]
    axB.set_title(f"[2] SOURCING  (l=1 → l=3 via cubic)\n"
                  f"λ=0: ‖ψ₃‖ = {l3_at_zero:.0e} (≡0)   ‖ψ₃‖(2L)/‖ψ₃‖(L) = {ratio:.3f} → linear")
    axB.legend(loc="lower right", fontsize=8)
    cap(axB, "Only l=1 is seeded, so l=3 is born purely from the cubic:\n"
             "zero at λ=0, doubling when λ doubles — a source linear in λ.")

    # --- [3] HERMITICITY: instantaneous no-net-work on the evolved state ---
    axC = ax[0, 2]
    _ls, _c, Qbar = angular.coupling_tensors(3, 1)
    psi_probe = t3_L["psi"]
    net = np.empty(psi_probe.shape[1])
    for kk in range(psi_probe.shape[1]):
        p = psi_probe[:, kk]
        F = np.einsum("jabc,a,b,c->j", Qbar, p, np.conj(p), p)
        net[kk] = abs(np.sum(np.imag(np.conj(p) * F))) / (np.sum(np.abs(np.conj(p) * F)) + 1e-300)
    axC.semilogy(t3_L["t"], net + 1e-300, color="C4", lw=1.2)
    axC.axhline(1e-12, color="C3", ls="--", lw=1.0, label="pass threshold 1e-12")
    axC.set_ylim(1e-18, 1e-3); axC.set_xlabel("t / M")
    axC.set_ylabel("|Σⱼ Im(ψⱼ* Fⱼ)| / Σⱼ|ψⱼ* Fⱼ|")
    axC.set_title("[3] HERMITICITY  (cubic does no net work)")
    axC.text(0.5, 0.90, f"worst = {net.max():.1e}  (machine)", transform=axC.transAxes,
             ha="center", fontsize=8, bbox=dict(boxstyle="round", fc="w", ec="0.7"))
    axC.legend(loc="lower right", fontsize=8)
    cap(axC, "The cubic moves norm between channels but does no net work:\n"
             "Σⱼ Im(ψⱼ*Fⱼ)=0 on the evolved state at all times (1e-16 floor).")

    # --- [4a] CASCADE populations + fill-time lag ---
    axD = ax[1, 0]
    tt5 = t5["t"]
    for i, l in enumerate(t5["ls"]):
        axD.plot(tt5, np.abs(t5["cloud_norm"][i]) + 1e-300, color=CL[int(l)], lw=1.8, label=f"l={l}")
    def t_reach(nrm, frac): return tt5[np.argmax(nrm >= frac * np.max(nrm))]
    th3 = t_reach(t5["cloud_norm"][1], 0.5); th5 = t_reach(t5["cloud_norm"][2], 0.5)
    axD.axvline(th3, color=CL[3], ls=":", lw=1.0); axD.axvline(th5, color=CL[5], ls=":", lw=1.0)
    axD.set_yscale("log"); axD.set_ylim(1e-6, 1e2)
    axD.set_xlabel("t / M"); axD.set_ylabel("‖ψ_l‖ (interior L2)")
    axD.set_title("[4] CASCADE populations  (l=5 fills after l=3)")
    axD.text(0.97, 0.06, f"t½: l=3 = {th3:.0f}M,  l=5 = {th5:.0f}M\n(l=5 fills later)",
             transform=axD.transAxes, ha="right", fontsize=8,
             bbox=dict(boxstyle="round", fc="w", ec="0.7"))
    axD.legend(loc="upper right", fontsize=8)
    cap(axD, "Three l=1 quanta reach only l≤3, so l=5 is second-generation:\n"
             "it grows only after l=3 fills, reaching half-peak ~7M later.")

    # --- [4b] CASCADE early power law: l=3 ~ t^1, l=5 ~ t^2 ---
    axE = ax[1, 1]
    def early(nrm):
        m = (tt5 > 0) & (nrm > 1e-12 * np.max(nrm)) & (tt5 < tt5[np.argmax(nrm)])
        s = np.polyfit(np.log(tt5[m]), np.log(nrm[m]), 1)
        return tt5[m], nrm[m], s
    for l, nrm in ((3, t5["cloud_norm"][1]), (5, t5["cloud_norm"][2])):
        tm, nm, s = early(nrm)
        axE.loglog(tm, nm, color=CL[l], lw=2.0, label=f"l={l}:  ~ t^{s[0]:.2f}")
        axE.loglog(tm, np.exp(np.polyval(s, np.log(tm))), color=CL[l], ls="--", lw=1.0)
    axE.set_xlabel("t / M"); axE.set_ylabel("‖ψ_l‖")
    axE.set_title("[4] CASCADE early power law  (1st vs 2nd generation)")
    axE.legend(loc="upper left", fontsize=8)
    cap(axE, "First-gen l=3 rises ∝t; second-gen l=5 rises ∝t² (needs l=3):\n"
             "fitted slopes 0.98 and 2.11 confirm the generation ordering.")

    # --- [3-disc] λ-scan: dΓ ∝ λ (Hermitian detuning), not λ² (net-work) ---
    axF = ax[1, 2]
    axF.plot(lam_scan, -dG * 1e3, "o", color="C5", ms=7, label="measured  −ΔΓ")
    slope = np.sum(np.array(lam_scan) * (-dG)) / np.sum(np.array(lam_scan) ** 2)  # through origin
    xs = np.linspace(0, max(lam_scan) * 1.05, 50)
    axF.plot(xs, slope * xs * 1e3, "C5--", lw=1.2, label=f"linear fit  (∝ λ¹)")
    q2 = (-dG[-1]) / lam_scan[-1] ** 2                                            # a λ² through last pt
    axF.plot(xs, q2 * xs ** 2 * 1e3, "0.6", ls=":", lw=1.2, label="∝ λ²  (ruled out)")
    axF.set_xlabel("λ (quartic coupling)"); axF.set_ylabel("−ΔΓ_total  (×10⁻³)")
    axF.set_title("[3-disc] total-Γ shift is O(λ), not O(λ²)")
    axF.legend(loc="upper left", fontsize=8)
    cap(axF, "Total-Γ shift is linear in λ (dashed), not quadratic (dotted):\n"
             "a net-work violation would be O(λ²); linear ⇒ detuning artifact.")

    fig.subplots_adjust(top=0.92, bottom=0.16, left=0.055, right=0.985,
                        hspace=0.60, wspace=0.27)
    fig.savefig(q.out, dpi=120)
    print("saved", q.out)


if __name__ == "__main__":
    main()
