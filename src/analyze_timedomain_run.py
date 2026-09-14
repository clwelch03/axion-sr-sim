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
from frequency_filter import fourier_amplitudes
from measure_gamma import measure_filtered_mode
import plotstyle  # sets the repo default colour cycle (C0..C9) on import

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


def plot_raw_spectrum(times, probe_signal, spin, alpha, m, output_path):
    """Save the raw probe power spectrum with the Leaver frequencies marked."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    angular_frequencies, amplitudes = fourier_amplitudes(times, probe_signal)
    # solver convention psi ~ exp(-i wR t) places the peak at -wR; flip to the
    # physical positive frequency for plotting.
    physical_frequency = -angular_frequencies
    power = np.abs(amplitudes) ** 2
    power /= power.max()

    colours = ["C0", "C1"]  # match fig13_growth_rates: n=0 -> C0, n=1 -> C1
    leaver_omega = {
        label: complex(leaver.solve(alpha, spin, l=m, m=m, n=n_p,
                                     spheroidal=True))
        for label, n_p in OVERTONES
    }

    figure, (wide_axis, zoom_axis) = plt.subplots(2, 1, figsize=(9, 8))
    for axis in (wide_axis, zoom_axis):
        axis.semilogy(physical_frequency, power, "k-", lw=0.7)
        for (label, _), colour in zip(OVERTONES, colours):
            axis.axvline(
                leaver_omega[label].real, color=colour, ls="--", lw=1.2,
                label=f"Leaver {label}: "
                      f"$\\omega M$ = {leaver_omega[label].real:.6f}")
        axis.set_ylabel(r"normalized raw power  $|f(\omega)|^2$")
        axis.legend(fontsize=8, loc="lower left")

    wide_axis.set_xlim(0.0, 0.5)
    wide_axis.set_ylim(1e-12, 2.0)
    wide_axis.set_title(
        f"Dolan Gaussian evolution: raw $\\ell=m={m}$ probe spectrum "
        f"($a={spin}M$, $M\\mu={alpha}$, $t_{{\\max}}={times[-1]:.0f}M$)")

    fundamental = leaver_omega["n=0"].real
    zoom_axis.set_xlim(fundamental - 0.014, fundamental + 0.016)
    zoom_axis.set_ylim(1e-8, 2.0)
    zoom_axis.set_xlabel(r"physical frequency  $\omega M$")

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    print(f"\nsaved {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("npz_path", help="run .npz written by measure_gamma.py")
    parser.add_argument("--search-half-width", type=float, default=0.002,
                        help="frequency window (in omega*M) around each Leaver "
                             "peak for the filtered-mode search")
    parser.add_argument("--spectrum-plot", default=None,
                        help="if given, save the raw probe power spectrum here")
    arguments = parser.parse_args()

    times, probe_signal, spin, alpha, m = load_probe_record(arguments.npz_path)
    print(f"loaded {times.size} samples, t_max={times[-1]:.0f}M, "
          f"a={spin}, alpha={alpha}, m={m}")
    report_growth_rates(times, probe_signal, spin, alpha, m,
                        arguments.search_half_width)
    if arguments.spectrum_plot is not None:
        plot_raw_spectrum(times, probe_signal, spin, alpha, m,
                          arguments.spectrum_plot)


if __name__ == "__main__":
    main()
