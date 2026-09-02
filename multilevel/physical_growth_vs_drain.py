#!/usr/bin/env python
"""Compare native superradiant growth against the native self-interaction
drain, using each mode's own Leaver Im(omega) as its linear rate (no hand
tuning), at the one alpha window where the intended sign pattern holds.

Step 0 (see module docstring notes below) found that at the demo's default
alpha=0.30, dimensionless_spin=0.99, MODE_322 is itself growing
(Im(omega_322) > 0), not damped, so it cannot act as a genuine linear sink.
Scanning dimensionless_spin=0.99 over alpha in [0.05, 0.80] shows
Im(omega_211) > 0 (211 superradiant) throughout alpha <~ 0.44, while
Im(omega_322) is negative (322 damped) only for alpha <~ 0.158 -- above that
322 itself turns superradiant. The overlap window is alpha in
[0.05, ~0.158]; we use alpha=0.15, comfortably inside it, as the working
point. At that point Im(omega_100) = -2.79e-4 is far larger in magnitude
than either Im(omega_211) ~ +1.8e-10 or Im(omega_322) ~ -4.3e-16 -- the
mediator's own decay, not the intended 322 sink, is the dominant native
linear rate in this three-mode system.

This uses ``transition_demo.run(..., use_leaver_linear_rates=True)``, the
opt-in extension added for this comparison; default behavior of
``transition_demo.run`` is unchanged (see that module's tests, still exactly
norm-conserving with the new flag left at its default ``False``).
"""
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from transition_demo import run as run_transition

ALPHA = 0.15
DIMENSIONLESS_SPIN = 0.99
FINAL_TIME = 3000.0
NATIVE_INTERACTION_STRENGTH = 1000.0  # transition_demo.run's own default
COMMON_KWARGS = dict(
    alpha=ALPHA,
    dimensionless_spin=DIMENSIONLESS_SPIN,
    final_time=FINAL_TIME,
    use_leaver_linear_rates=True,
)


def log_amplitude_slope(t, amplitude_211):
    """Least-squares slope of ln|A_211| vs t: the net secular growth/decay
    rate, averaged over any Rabi-like oscillation from the cubic coupling."""
    log_magnitude = np.log(np.abs(amplitude_211))
    slope, _intercept = np.polyfit(t, log_magnitude, 1)
    return slope


def main():
    print(f"alpha={ALPHA}, dimensionless_spin={DIMENSIONLESS_SPIN}, "
          f"final_time={FINAL_TIME}")

    start_time = time.time()
    pure_linear_control = run_transition(
        interaction_strength=0.0, **COMMON_KWARGS
    )
    print(f"pure-linear control run: {time.time() - start_time:.1f} s")

    start_time = time.time()
    native_coupling_run = run_transition(
        interaction_strength=NATIVE_INTERACTION_STRENGTH, **COMMON_KWARGS
    )
    print(f"native-coupling run: {time.time() - start_time:.1f} s")

    labels = [mode.label for mode in pure_linear_control["active_modes"]]
    index_211 = labels.index("211")
    index_322 = labels.index("322")
    index_100 = labels.index("100")
    complex_frequencies = pure_linear_control["complex_frequencies"]
    imag_211 = complex_frequencies[index_211].imag
    imag_322 = complex_frequencies[index_322].imag
    imag_100 = complex_frequencies[index_100].imag
    print("active modes:", labels)
    print(f"Im(omega_211) = {imag_211:.6e}  (native superradiant growth rate)")
    print(f"Im(omega_322) = {imag_322:.6e}  (native linear sink)")
    print(f"Im(omega_100) = {imag_100:.6e}  (mediator's own decay)")

    t = pure_linear_control["t"]
    control_amplitude_211 = pure_linear_control["amplitudes"][:, index_211]
    native_amplitude_211 = native_coupling_run["amplitudes"][:, index_211]

    control_rate = log_amplitude_slope(t, control_amplitude_211)
    net_rate = log_amplitude_slope(t, native_amplitude_211)
    print(f"pure-linear control: fitted d ln|A_211|/dt = {control_rate:.4e} "
          f"(expected Im(omega_211) = {imag_211:.4e})")
    print(f"native-coupling run: fitted net d ln|A_211|/dt = {net_rate:.4e}")
    print(f"ratio net_rate / Im(omega_211) = {net_rate / imag_211:.3f}")

    control_populations = pure_linear_control["populations"]
    native_populations = native_coupling_run["populations"]
    print("native run max populations:",
          dict(zip(labels, native_populations.max(axis=0))))
    print("native run final populations:",
          dict(zip(labels, native_populations[-1])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"211": "C0", "100": "C1", "322": "C2"}

    # Panel A: fractional |A_211|^2 deviation from its initial value, pure
    # linear growth only vs linear growth plus native self-interaction.
    initial_population_211 = control_populations[0, index_211]
    control_fractional_deviation = (
        control_populations[:, index_211] / initial_population_211 - 1.0
    )
    native_fractional_deviation = (
        native_populations[:, index_211] / initial_population_211 - 1.0
    )
    ax[0].plot(t, control_fractional_deviation, color=colors["211"], lw=2,
               ls="--", label=r"linear only ($Im(\omega_{211})$, $\lambda=0$)")
    ax[0].plot(t, native_fractional_deviation, color=colors["211"], lw=1.5,
               label=r"linear + native self-interaction ($\lambda=1000$)")
    ax[0].axhline(0.0, color="0.7", lw=0.8)
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"$|A_{211}|^2 / |A_{211}(0)|^2 - 1$")
    ax[0].set_title(r"$|211\rangle$ population: native growth vs native drain")
    ax[0].legend(fontsize=9)

    # Panel B: native linear rates (bars) vs the extracted net |211> rate
    # under the coupling, all on a common log|rate| scale with sign shown.
    rate_names = [
        r"$Im(\omega_{211})$" + "\n(211 growth)",
        r"$Im(\omega_{322})$" + "\n(322 sink)",
        r"$Im(\omega_{100})$" + "\n(mediator decay)",
        r"net $d\ln|A_{211}|/dt$" + "\n(with coupling)",
    ]
    rate_values = [imag_211, imag_322, imag_100, net_rate]
    bar_colors = ["C0" if value >= 0 else "C3" for value in rate_values]
    magnitudes = np.abs(rate_values)
    bars = ax[1].bar(rate_names, magnitudes, color=bar_colors)
    ax[1].set_yscale("log")
    ax[1].set_ylabel(r"$|{\rm rate}|$  (blue = growth, red = decay)")
    ax[1].set_title("native rates at "
                     rf"$\alpha={ALPHA}$, $\chi={DIMENSIONLESS_SPIN}$")
    for bar, value in zip(bars, rate_values):
        ax[1].annotate(f"{value:+.2e}", xy=(bar.get_x() + bar.get_width() / 2,
                        bar.get_height()), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8)

    fig.suptitle(
        "physically calibrated linear rates (no hand tuning): native\n"
        "self-interaction does not measurably compete with native growth"
        f"\n(net rate / bare 211 growth rate = {net_rate / imag_211:.2f}x, "
        "both far below any accessible timescale)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    output_path = os.path.join(_HERE, "physical_growth_vs_drain.png")
    fig.savefig(output_path, dpi=130)
    print("saved", output_path)


if __name__ == "__main__":
    main()
