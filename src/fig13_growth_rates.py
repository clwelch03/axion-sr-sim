"""Reproduce Dolan 1212.1477 Fig. 13: superradiant growth rate M*nu vs mass
coupling alpha at a=0.99M, for the l=m=1 fundamental (n=0) and first overtone
(n=1).

THREE THINGS ON ONE SET OF AXES, and the difference between them is the point:

1. Leaver (dashed lines) -- leaver.py's Dolan-2007 continued fraction, the
   independent frequency-domain benchmark both of our methods are measured
   against.
2. Our eigenvalue solve (open circles) -- solve_qbs diagonalises the discretised
   linear generator and reads the QBS complex frequency straight off. Cheap and
   accurate, but it DOES NOT EXTEND TO THE NONLINEAR PROBLEM: once the field
   self-interacts (|psi|^2 psi) or backreacts there is no linear generator to
   diagonalise.
3. Our time-domain measurement (filled squares) -- Dolan's Gaussian initial data
   evolved to t=10^6 M, then the frequency filter of Sec. IV A 2 (Eq. 46/49)
   applied to the fixed-x probe series. No linear generator appears anywhere in
   that chain, so this is the path that carries over to the nonlinear problem.
   Each point is one heavy ``measure_gamma.py`` run; see ``runs/README.md``.

So (2) validates the discretisation and (3) validates the extraction method we
actually intend to use later. Agreement of (3) with (1) is what licenses
trusting the time-domain rate when there is no Leaver curve to check against.

Overtone convention: leaver.py's principal number n_p = l+1+n_radial, so the
Fig.-13 fundamental n=0 is leaver n_p=2 and the overtone n=1 is n_p=3.

    python src/fig13_growth_rates.py                     # both methods
    python src/fig13_growth_rates.py --skip-eigenvalue   # Leaver + time domain
"""
import argparse
import os
import time

import mpmath as mp
import numpy as np
import leaver
from scipy.sparse.linalg import ArpackNoConvergence
from solver import Tower
from qbs_eigenmode import solve_qbs
import fig13_timedomain_points
import plotstyle  # sets the repo default colour cycle (C0..C9) on import

SPIN = 0.99
curve_alpha_values = np.linspace(0.30, 0.44, 71)  # d_alpha=0.002
point_alpha_values = np.round(np.arange(0.30, 0.441, 0.01), 3)
# (Fig.-13 label, leaver principal number, plot colour); C0/C1 are the first
# two repo brand colours (see plotstyle.py).
overtones = [("n=0", 2, "C0"), ("n=1", 3, "C1")]
ARPACK_MAXITER = 2000
MAX_FIG13_RATE = 1.0e-4


def our_eigenvalue(alpha, leaver_omega):
    """Growth rate from our discretised-operator eigenvalue solve, targeting the
    full complex Leaver frequency so the right overtone is selected."""
    tower = Tower(spin=SPIN, alpha=alpha, m=1, l_max=1,
                  x_left=-200.0, x_right=250.0, grid_spacing_n=8,
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


def leaver_curves():
    """Track both Fig.-13 branches across alpha.  Returns {label: omega array}."""
    curves = {}
    for label, n_p, _ in overtones:
        print(f"\n--- {label} (leaver n_p={n_p}), a={SPIN}: "
              f"tracking Leaver curve ---", flush=True)
        curves[label] = leaver_curve(n_p)
    return curves


def leaver_omega_by_alpha(curves=None):
    """Branch-tracked Leaver frequency at each point alpha, {label: {alpha: w}}.

    This, not a per-alpha ``leaver.solve``, is the right reference for the
    time-domain extraction: an unguided continued-fraction solve leaves the
    physical branch near the top of the sweep (at alpha=0.44 the n=1 root finder
    returns Im(omega) ~ -0.1, which is no quasi-bound state).
    """
    if curves is None:
        curves = leaver_curves()
    return {
        label: {round(alpha, 3): omega
                for alpha, omega in zip(curve_alpha_values, curves[label])}
        for label in curves
    }


def compute(curves=None):
    """Leaver curves plus our eigenvalue points at each Fig.-13 point alpha."""
    if curves is None:
        curves = leaver_curves()
    omega_by_alpha = leaver_omega_by_alpha(curves)
    results = {}
    for label, n_p, _ in overtones:
        curve_omega = curves[label]
        curve_by_alpha = omega_by_alpha[label]
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


# Distinct marker per grid resolution, so several dx overlay on one figure
# without ambiguity.  Colour still encodes the overtone (red n=0, blue n=1).
GRID_MARKERS = {4: "s", 6: "D", 8: "^", 12: "v"}


def plot(curves, results, timedomain_sets, path):
    """Draw Fig. 13: Leaver curves, eigenvalue points, time-domain points.

    ``results`` is optional (None to omit the eigenvalue points).
    ``timedomain_sets`` is a list of ``(grid_spacing_n, (alpha_values,
    measured))`` -- one entry per time-domain resolution to overlay.  Each grid
    gets its own marker and a legend label naming dx = M / grid_spacing_n.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5.25))  # 4:3
    for label, n_p, colour in overtones:
        ax.plot(curve_alpha_values, curves[label].imag, "--", color=colour,
                lw=1.3, label=f"{label}  Leaver (continued fraction)")
        if results is not None:
            _, our_omega = results[label]
            ax.plot(point_alpha_values, our_omega.imag, "o", color=colour, ms=6,
                    mfc="none", label=f"{label}  our eigenvalue solve")
        for grid_spacing_n, (alpha_values, measured) in timedomain_sets:
            _, measured_gamma, _, _ = measured[label]
            marker = GRID_MARKERS.get(grid_spacing_n, "s")
            ax.plot(alpha_values, measured_gamma, marker, color=colour, ms=5,
                    label=fr"{label}  time domain, $\Delta x = M/{grid_spacing_n}$")

    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_ylim(bottom=-0.1e-7)  # crop the deep n=4 n=1 outlier at alpha=0.44
    ax.set_xlabel(r"mass coupling  $M\mu = \alpha$")
    ax.set_ylabel(r"growth rate  $M\nu = \mathrm{Im}\,\omega$")
    ax.set_title(f"Superradiant growth rate on Kerr, a={SPIN}M  (l=m=1)")
    ax.legend(fontsize=8, loc="upper left")
    # ax.text(0.98, 0.02,
    #         "open circles: eigenvalue solve -- validates the LINEAR rates,\n"
    #         "but has no nonlinear extension (no generator to diagonalise)\n"
    #         "filled squares: time domain + frequency filter -- slower, but\n"
    #         "this is the method that carries over to the nonlinear problem",
    #         transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
    #         color="0.35", style="italic")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"\nsaved {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-eigenvalue", action="store_true",
                        help="omit the eigenvalue points, i.e. skip the slow "
                             "ARPACK sweep and redraw from the Leaver curves "
                             "and the finished time-domain runs alone")
    parser.add_argument("--grid-spacing-n", type=int, nargs="+", default=[4],
                        help="grid resolution(s) of the time-domain runs to "
                             "overlay; pass several to compare, e.g. "
                             "--grid-spacing-n 4 8")
    arguments = parser.parse_args()

    curves = leaver_curves()
    results = None if arguments.skip_eigenvalue else compute(curves)

    omega_by_alpha = leaver_omega_by_alpha(curves)
    timedomain_sets = []
    for grid_spacing_n in arguments.grid_spacing_n:
        run_records = fig13_timedomain_points.find_run_records(
            fig13_timedomain_points.DEFAULT_RUN_DIRECTORY, grid_spacing_n)
        if not run_records:
            print(f"\n--- time domain: no grid_spacing_n={grid_spacing_n} runs, "
                  f"skipping ---", flush=True)
            continue
        print(f"\n--- time domain: {len(run_records)} runs at "
              f"grid_spacing_n={grid_spacing_n} ---", flush=True)
        alpha_values, measured = fig13_timedomain_points.growth_rate_points(
            run_records, omega_by_alpha)
        fig13_timedomain_points.print_growth_rate_table(alpha_values, measured)
        timedomain_sets.append((grid_spacing_n, (alpha_values, measured)))

    output_path = "plots/fig13_growth_rates.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plot(curves, results, timedomain_sets, output_path)


if __name__ == "__main__":
    main()
