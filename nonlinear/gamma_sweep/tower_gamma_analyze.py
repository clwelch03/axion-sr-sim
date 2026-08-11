#!/usr/bin/env python
"""Analyze the Gamma(A, l_max) growth-suppression sweep (tower_gamma_sweep.sbatch).

Globs the tower_gamma_lmax*_A*_lam*.npz outputs, re-extracts the l=1 growth rate
with freqfilter (offline, so kbins is tunable), differences against the lambda=0
baseline to cancel the freqfilter/PML systematic, and plots Gamma_l1(A) per l_max.

Prediction under test: l_max=1 is the CONTROL (Gamma flat at +1.5e-7, Hermitian
cubic => detuning only); l_max>=3 FALLS as the conservative cubic drains the
growing l=1 into the strongly-damped l=3,5 -> horizon. Runs fine on a PARTIAL
sweep (skips files not yet written), so it can be run while jobs are still going.

    gwenv/bin/python tower_gamma_analyze.py [--kbins 5] [--plot tower_gamma_results.png]
"""
import sys, os, glob, argparse, numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_RD = os.path.dirname(os.path.dirname(_HERE))            # radial1d/
sys.path.insert(0, _RD)
import evolve as ev


def spectral_peak(t, psi):
    """Dominant frequency E of psi ~ e^{-iEt} (np.fft peak sits at omega=-E).
    Used to RE-CENTER freqfilter on the amplitude-DETUNED peak: at finite A the
    Hermitian cubic shifts E up toward m*Omega_H, and a filter pinned to the
    linear Leaver E clips the shifted mode and returns garbage (seen at A>=2)."""
    N = len(t); w = np.fft.fft(psi)
    omega = 2 * np.pi * np.fft.fftfreq(N, d=(t[-1] - t[0]) / (N - 1))
    return -omega[np.argmax(np.abs(w))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=_HERE, help="directory of tower_gamma_*.npz")
    ap.add_argument("--kbins", type=float, default=5.0, help="freqfilter half-width (Dolan ~5)")
    ap.add_argument("--plot", default=None)
    q = ap.parse_args()

    files = sorted(glob.glob(os.path.join(q.dir, "tower_gamma_lmax*_A*_lam*.npz")))
    if not files:
        print(f"no tower_gamma_*.npz in {q.dir} yet (jobs still running?)"); return

    rows = []          # (l_max, A, lam, E_l1, G_l1, per-channel dict)
    GL = None
    for f in files:
        d = np.load(f, allow_pickle=True)
        w = complex(d["w_leaver"]); GL = w.imag
        t = d["t"]; psi = d["psi"]; ls = list(np.atleast_1d(d["ls"]))
        chan = {}
        for i, l in enumerate(ls):
            E, G = ev.freqfilter(t, psi[i], E_peak=spectral_peak(t, psi[i]), kbins=q.kbins)
            chan[int(l)] = (E, G)
        E1, G1 = chan[1]
        rows.append(dict(l_max=int(d["l_max"]), A=float(d["amp"]), lam=float(d["lam"]),
                         E1=E1, G1=G1, chan=chan, nfile=os.path.basename(f)))

    # baseline: the lambda=0 run's l=1 rate (freqfilter systematic anchor)
    base = [r for r in rows if r["lam"] == 0.0]
    G_base = base[0]["G1"] if base else GL
    print(f"Leaver Gamma_l1 = {GL:+.4e}   freqfilter baseline (lambda=0) = {G_base:+.4e} "
          f"({G_base/GL:+.3f}x)   [{len(rows)}/{10} runs present, kbins={q.kbins}]\n")

    print(f"{'l_max':>6}{'A':>6}{'lambda':>8}{'E_l1':>11}{'Gamma_l1':>13}"
          f"{'ratio':>9}{'dGamma=G-Gbase':>16}", flush=True)
    for r in sorted(rows, key=lambda r: (r["l_max"], r["lam"], r["A"])):
        dG = r["G1"] - G_base
        print(f"{r['l_max']:>6}{r['A']:>6g}{r['lam']:>8g}{r['E1']:>11.6f}"
              f"{r['G1']:>13.4e}{r['G1']/GL:>8.3f}x{dG:>16.4e}")

    if q.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        nl = sorted({r["l_max"] for r in rows if r["lam"] != 0.0})
        fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
        CL = {1: "C0", 3: "C1", 5: "C2", 7: "C3"}
        for lm in nl:
            pts = sorted([r for r in rows if r["l_max"] == lm and r["lam"] != 0.0],
                         key=lambda r: r["A"])
            if not pts: continue
            A = [r["A"] for r in pts]; G = [r["G1"] for r in pts]
            ax[0].plot(A, G, "o-", color=CL.get(lm, "k"), label=f"l_max={lm}")
            dG = [abs(r["G1"] - G_base) for r in pts]
            ax[1].plot(A, dG, "o-", color=CL.get(lm, "k"), label=f"l_max={lm}")
        ax[0].axhline(GL, color="k", ls="--", lw=1, label=f"Leaver +{GL:.2e}")
        ax[0].axhline(0, color="0.5", lw=0.8)
        ax[0].set_xlabel("seed amplitude A"); ax[0].set_ylabel("Gamma_l1 (freqfilter)")
        ax[0].set_title("l=1 growth rate vs amplitude\n(l_max curves coincide => detuning-dominated, not multi-mode)")
        ax[0].legend(fontsize=8)
        ax[1].set_xscale("log"); ax[1].set_yscale("log")
        ax[1].set_xlabel("seed amplitude A"); ax[1].set_ylabel("|ΔΓ| = |Γ(A) − Γ_baseline|")
        ax[1].set_title("growth suppression vs A\n(slope 2 ⇒ ΔΓ ∝ A²; l_max convergence)")
        ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(q.plot, dpi=120); print("\nsaved", q.plot)


if __name__ == "__main__":
    main()
