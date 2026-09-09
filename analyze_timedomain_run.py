"""Extract growth rates from a finished Dolan time-domain run.

Loads a probe record written by ``measure_gamma.py`` (the fixed-x complex time
series of the lowest angular channel) and reports, for each Fig.-13 overtone,
the frequency-filter growth rate Gamma (Dolan Eq. 49) against the runtime
Leaver benchmark. Captures, as a rerunnable script, the Gamma table quoted in
the README.

    python analyze_timedomain_run.py runs/dolan_fig13_a099_alpha042_t100000.npz
"""

import argparse

import numpy as np

import leaver
from measure_gamma import measure_filtered_mode

# Dolan Fig. 13 tracks the l=m=1 fundamental (n=0) and first overtone (n=1).
# leaver.py's principal number is n_p = l + 1 + n_radial, so n=0 -> n_p=2.
OVERTONES = [("n=0", 2), ("n=1", 3)]


def load_probe_record(npz_path):
    """Return (times, probe_signal, spin, alpha, m) from a run .npz."""
    with np.load(npz_path) as record:
        times = record["sample_times"]
        probe_signal = record["probe_signal"]
        spin = float(record["tower_parameters"][0])
        alpha = float(record["tower_parameters"][1])
        m = int(record["tower_parameters"][2])
    return times, probe_signal, spin, alpha, m


def report_growth_rates(times, probe_signal, spin, alpha, m,
                        search_half_width=0.002):
    """Print measured (filtered) Gamma vs spheroidal Leaver for each overtone."""
    print(f"\n{'mode':5s} {'wR (Leaver)':>13s} {'wR (meas)':>13s}"
          f" {'G (Leaver)':>13s} {'G (meas)':>13s} {'ratio':>7s}")
    for label, n_p in OVERTONES:
        reference = complex(leaver.solve(alpha, spin, l=m, m=m, n=n_p,
                                         spheroidal=True))
        measured_frequency, measured_gamma, _ = measure_filtered_mode(
            times, probe_signal, reference.real, search_half_width)
        ratio = measured_gamma / reference.imag if reference.imag else np.nan
        print(f"{label:5s} {reference.real:13.6f} {measured_frequency:13.6f}"
              f" {reference.imag:13.4e} {measured_gamma:13.4e} {ratio:7.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("npz_path", help="run .npz written by measure_gamma.py")
    parser.add_argument("--search-half-width", type=float, default=0.002,
                        help="frequency window (in omega*M) around each Leaver "
                             "peak for the filtered-mode search")
    arguments = parser.parse_args()

    times, probe_signal, spin, alpha, m = load_probe_record(arguments.npz_path)
    print(f"loaded {times.size} samples, t_max={times[-1]:.0f}M, "
          f"a={spin}, alpha={alpha}, m={m}")
    report_growth_rates(times, probe_signal, spin, alpha, m,
                        arguments.search_half_width)


if __name__ == "__main__":
    main()
