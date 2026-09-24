"""Time-domain Fig. 13 points: the frequency-filter Gamma of each finished run.

This is the OTHER half of the Fig. 13 reproduction.  ``fig13_growth_rates.py``
gets the growth rate by diagonalising the linear generator; that is cheap and
accurate but stops at the linear problem.  Here the same rate comes out of a
full time-domain evolution -- Dolan's Gaussian initial data (Eq. 45) evolved to
t=10^6 M, then the frequency filter (Eq. 46/49) applied to the fixed-x probe
series.  That path has no linear generator in it, so it is the one that carries
over to the nonlinear problem.

Each alpha is one ``measure_gamma.py`` run under ``runs/dolan/fig13/sweep/``; the array
recipe is in ``runs/README.md``.  A run stores its own alpha and grid
resolution in ``tower_parameters``, so runs are selected by reading the files,
not by parsing their names.

The Leaver reference frequency is supplied by the caller rather than solved per
alpha, because a per-alpha continued-fraction solve leaves the physical branch
near the top of the sweep -- at alpha=0.44 the n=1 root finder lands on
Im(omega) ~ -0.1, which is not a quasi-bound state at all.  The branch-tracked
curve in ``fig13_growth_rates.leaver_curve`` is the reliable source.

    PYTHONPATH=src python -m dolan.fig13_timedomain_points --grid-spacing-n 4
"""

import argparse
import os

import numpy as np

from dolan.analyze_timedomain_run import load_probe_record
from dolan.measure_gamma import measure_filtered_mode

# Dolan Fig. 13 tracks the l=m=1 fundamental (n=0) and first overtone (n=1).
# leaver.py's principal number is n_p = l + 1 + n_radial, so n=0 -> n_p=2.
OVERTONES = [("n=0", 2), ("n=1", 3)]
DEFAULT_RUN_DIRECTORY = "runs/dolan/fig13/sweep"
DEFAULT_SEARCH_HALF_WIDTH = 0.002


def find_run_records(run_directory, grid_spacing_n):
    """Return (alpha, path) for every run in ``run_directory`` at this grid n.

    Selection reads ``tower_parameters`` out of each .npz -- the run's own
    record of what it was -- so a renamed or re-tagged file cannot be
    misfiled.  Results are sorted by alpha.
    """
    records = []
    for entry in sorted(os.listdir(run_directory)):
        if not entry.endswith(".npz"):
            continue
        path = os.path.join(run_directory, entry)
        with np.load(path) as record:
            tower_parameters = record["tower_parameters"]
        if int(tower_parameters[3]) != grid_spacing_n:
            continue
        records.append((float(tower_parameters[1]), path))
    return sorted(records)


def growth_rate_points(run_records, leaver_omega_by_alpha,
                       search_half_width=DEFAULT_SEARCH_HALF_WIDTH):
    """Measure the filtered Gamma of every run, for every Fig.-13 overtone.

    ``leaver_omega_by_alpha`` maps ``label -> {rounded alpha -> complex omega}``
    and supplies both the search centre (its real part) and the benchmark (its
    imaginary part).  A run whose alpha is missing from that map, or whose
    filter search fails, records NaN rather than aborting the sweep: one bad
    alpha should not cost the other fourteen.

    Returns ``alpha_values`` and a dict ``label -> (leaver_gamma,
    measured_gamma, measured_frequency, final_time)``, all NumPy arrays.
    """
    alpha_values = np.array([alpha for alpha, _ in run_records])
    results = {}
    for label, _ in OVERTONES:
        leaver_gamma = np.full(alpha_values.size, np.nan)
        measured_gamma = np.full(alpha_values.size, np.nan)
        measured_frequency = np.full(alpha_values.size, np.nan)
        final_time = np.full(alpha_values.size, np.nan)

        for index, (alpha, path) in enumerate(run_records):
            reference = leaver_omega_by_alpha[label].get(round(alpha, 3))
            if reference is None:
                continue
            leaver_gamma[index] = reference.imag

            times, probe_signal, _, _, _ = load_probe_record(path)
            final_time[index] = times[-1]
            try:
                frequency, gamma, _ = measure_filtered_mode(
                    times, probe_signal, reference.real, search_half_width)
            except ValueError:
                continue
            measured_frequency[index] = frequency
            measured_gamma[index] = gamma

        results[label] = (leaver_gamma, measured_gamma,
                          measured_frequency, final_time)
    return alpha_values, results


def print_growth_rate_table(alpha_values, results):
    """Print measured vs branch-tracked Leaver Gamma, one row per alpha."""
    print(f"\n{'alpha':>6s} {'t_final':>9s}", end="")
    for label, _ in OVERTONES:
        print(f" | {label + ' leaver':>13s} {label + ' meas':>13s} {'ratio':>6s}",
              end="")
    print()
    for index, alpha in enumerate(alpha_values):
        final_time = results[OVERTONES[0][0]][3][index]
        print(f"{alpha:6.2f} {final_time:9.0f}", end="")
        for label, _ in OVERTONES:
            leaver_gamma, measured_gamma, _, _ = results[label]
            ratio = (measured_gamma[index] / leaver_gamma[index]
                     if leaver_gamma[index] else np.nan)
            print(f" | {leaver_gamma[index]:13.4e} {measured_gamma[index]:13.4e}"
                  f" {ratio:6.3f}", end="")
        print(flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-directory", default=DEFAULT_RUN_DIRECTORY,
                        help="directory of measure_gamma.py run .npz files")
    parser.add_argument("--grid-spacing-n", type=int, default=4,
                        help="select only runs at this grid resolution")
    parser.add_argument("--search-half-width", type=float,
                        default=DEFAULT_SEARCH_HALF_WIDTH,
                        help="frequency window (in omega*M) around each Leaver "
                             "peak for the filtered-mode search")
    arguments = parser.parse_args()

    # Imported here, not at module scope: fig13_growth_rates imports this
    # module for the overlay, so a top-level import would be circular.
    from dolan.fig13_growth_rates import leaver_omega_by_alpha

    run_records = find_run_records(arguments.run_directory,
                                   arguments.grid_spacing_n)
    if not run_records:
        raise SystemExit(
            f"no grid_spacing_n={arguments.grid_spacing_n} runs in "
            f"{arguments.run_directory}")
    print(f"{len(run_records)} runs at grid_spacing_n="
          f"{arguments.grid_spacing_n} in {arguments.run_directory}",
          flush=True)

    alpha_values, results = growth_rate_points(
        run_records, leaver_omega_by_alpha(), arguments.search_half_width)
    print_growth_rate_table(alpha_values, results)


if __name__ == "__main__":
    main()
