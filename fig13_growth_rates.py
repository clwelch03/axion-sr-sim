"""Reproduce Dolan 1212.1477 Fig. 13: superradiant growth rate M*nu vs mass
coupling alpha at a=0.99M, for the l=m=1 fundamental (n=0) and first overtone
(n=1).

METHOD (read this before trusting the numbers): the points here come from our
1+1D *eigenvalue* solve -- solve_qbs diagonalises the discretised linear
generator and reads off the QBS complex frequency directly. That is a
frequency-domain method. It is cheap and accurate for the LINEAR problem, but it
DOES NOT EXTEND TO THE NONLINEAR PROBLEM: once the field self-interacts
(|psi|^2 psi) or backreacts, there is no linear generator to diagonalise, and
the growth rate must instead be extracted from a time-domain evolution
(DMD / Dolan's frequency-filter method, Sec. IV A 2 -- to be built later).
So this figure validates the linear rates and our machinery against Leaver; it
is not the tool for the eventual nonlinear science.

Lines: leaver.py (Dolan-2007 continued fraction), the independent frequency-domain
benchmark. Overtone convention: leaver.py's principal number n_p = l+1+n_radial,
so the Fig.-13 fundamental n=0 is leaver n_p=2 and the overtone n=1 is n_p=3.
"""
import time

import mpmath as mp
import numpy as np
import leaver
from scipy.sparse.linalg import ArpackNoConvergence
from solver import Tower
from qbs_eigenmode import solve_qbs

SPIN = 0.99
curve_alpha_values = np.linspace(0.30, 0.44, 71)  # d_alpha=0.002
point_alpha_values = np.round(np.arange(0.30, 0.441, 0.01), 3)
# (Fig.-13 label, leaver principal number, plot colour)
overtones = [("n=0", 2, "tab:red"), ("n=1", 3, "tab:blue")]
ARPACK_MAXITER = 2000
MAX_FIG13_RATE = 1.0e-4


def our_eigenvalue(alpha, leaver_omega):
    """Growth rate from our discretised-operator eigenvalue solve, targeting the
    full complex Leaver frequency so the right overtone is selected."""
    tower = Tower(spin=SPIN, alpha=alpha, m=1, l_max=1,
                  x_left=-200.0, x_right=250.0, n=8,
                  x_pml_center=-100.0, pml_width=10.0)
    omega, _ = solve_qbs(tower, omega_guess=leaver_omega,
                         maxiter=ARPACK_MAXITER)
    return omega


def leaver_curve(principal_number):
    """Track one QBS branch across alpha with a two-point linear predictor."""
    roots = []
    for index, alpha in enumerate(curve_alpha_values):
        if index < 2:
            root = leaver.solve(alpha, SPIN, l=1, m=1,
                                n=principal_number, spheroidal=True)
        else:
            alpha_step = curve_alpha_values[index - 1] - curve_alpha_values[index - 2]
            root_slope = (roots[-1] - roots[-2]) / alpha_step
            predicted_root = roots[-1] + root_slope * (
                alpha - curve_alpha_values[index - 1])
            root = leaver.solve(
                alpha, SPIN, l=1, m=1, n=principal_number, spheroidal=True,
                guess=mp.mpc(predicted_root.real, predicted_root.imag))

        if not 0.0 < root.real < alpha or abs(root.imag) > MAX_FIG13_RATE:
            raise RuntimeError(
                f"Leaver left the physical Fig. 13 branch at alpha={alpha:.3f}: "
                f"omega={root}")
        roots.append(root)
    return np.asarray(roots)


def compute():
    results = {}
    for label, n_p, _ in overtones:
        print(f"\n--- {label} (leaver n_p={n_p}), a={SPIN}: "
              f"tracking Leaver curve ---", flush=True)
        curve_omega = leaver_curve(n_p)
        curve_by_alpha = {
            round(alpha, 3): omega
            for alpha, omega in zip(curve_alpha_values, curve_omega)
        }
        our_omega = []
        print(f"\n--- {label} (leaver n_p={n_p}), a={SPIN} ---", flush=True)
        for alpha in point_alpha_values:
            point_start = time.perf_counter()
            w_l = curve_by_alpha[round(alpha, 3)]

            try:
                w_o = our_eigenvalue(alpha, w_l)
            except ArpackNoConvergence:
                elapsed_seconds = time.perf_counter() - point_start
                w_o = complex(np.nan, np.nan)
                print(f"  alpha={alpha:.2f}  ARPACK did not converge "
                      f"({elapsed_seconds:.1f} s); recording NaN", flush=True)
                our_omega.append(w_o)
                continue

            elapsed_seconds = time.perf_counter() - point_start
            our_omega.append(w_o)
            ratio = w_o.imag / w_l.imag if w_l.imag != 0 else float("nan")
            print(f"  alpha={alpha:.2f}  leaver nu={w_l.imag:+.4e}  "
                  f"ours nu={w_o.imag:+.4e}  ({ratio:.4f}x)  "
                  f"[{elapsed_seconds:.1f} s]", flush=True)
        results[label] = (curve_omega, np.asarray(our_omega))
    return results


def plot(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, n_p, colour in overtones:
        curve_omega, our_omega = results[label]
        ax.plot(curve_alpha_values, curve_omega.imag, "--", color=colour, lw=1.3,
                label=f"{label}  Leaver (continued fraction)")
        ax.plot(point_alpha_values, our_omega.imag, "o", color=colour, ms=5,
                mfc="none", label=f"{label}  our eigenvalue solve")

    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel(r"mass coupling  $M\mu = \alpha$")
    ax.set_ylabel(r"growth rate  $M\nu = \mathrm{Im}\,\omega$")
    ax.set_title(f"Superradiant growth rate on Kerr, a={SPIN}M  (l=m=1)")
    ax.legend(fontsize=8, loc="upper left")
    # the point the user asked to emphasise, on the figure itself
    ax.text(0.98, 0.02,
            "frequency-domain (eigenvalue) method:\nvalidates the LINEAR rates; "
            "does NOT extend to the\nnonlinear problem (no generator to "
            "diagonalise there)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color="0.35", style="italic")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"\nsaved {path}")


if __name__ == "__main__":
    results = compute()
    plot(results, "plots/fig13_growth_rates.png")
