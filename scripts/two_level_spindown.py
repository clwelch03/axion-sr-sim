"""Perturbative two-level spin-down: the rate-equation model of Baryakhtar, Galanis, Lasenby & Simon
(arXiv:2011.11646, docs/papers/2011.11646.pdf) for the |211> and |322> levels, as a comparison for the nonlinear
evolution of scripts/evolution_spindown.jl.

The model treats each level as a hydrogenic cloud with a number of quanta, and the quartic self-interaction as two
scattering processes at leading order (their Eqs. 34, 35 and 37, the terms outside the parentheses, with the rate
constants of their Table I):
  - 211 x 211 -> 322 x BH: two |211> quanta become one |322> quantum and one quantum absorbed by the hole,
        rate constant 4e-7 alpha^11 (M_Pl/f)^4 (1 + sqrt(1 - a^2)) mu;
  - 322 x 322 -> 211 x infinity: two |322> quanta become one |211> quantum and one quantum radiated to infinity,
        rate constant 1e-8 alpha^8 (M_Pl/f)^4 mu.
Both are second order in the quartic coupling. Their coupling is that of the cosine potential; for a potential whose
quartic term is `quartic_ratio` times the cosine's, both rates are multiplied by quartic_ratio^2. For the chiral
potential quartic_ratio = 1 - 3 z / (1 + z)^2 with z = m_u/m_d (0.3426 for the code's default z = 0.48).

Differences from the paper's equations, all on the linear side: the superradiant rates and the real frequencies are
the exact linear (spheroidal Leaver) ones, not the leading small-alpha expressions, and the hole's mass and angular
momentum are evolved with the energy and angular momentum of every quantum the hole exchanges.

Units and state as in scripts/leaver_spindown.py: time in M0, occupations N_i in f_a^2 M0^2 (script variables
`occupations`), b = M/M0, j = J/M0^2, alpha = alpha0 b, a = j/b^2, kappa = (f_a/M_Pl)^2. The paper's normalized
occupation is epsilon_i = N_i^physical / (G M^2) = kappa N_i / b^2. In these units the two scattering terms do not
contain f_a:

    R_BH       = quartic_ratio^2 * 4e-7 * alpha^11 * (1 + sqrt(1 - a^2)) * alpha0 * N_211^2 N_322 / b^4
    R_infinity = quartic_ratio^2 * 1e-8 * alpha^8 * alpha0 * N_322^2 N_211 / b^4
    dN_211/dt = 2 Im(omega_211) N_211 - 2 R_BH + R_infinity
    dN_322/dt = 2 Im(omega_322) N_322 + R_BH - 2 R_infinity
    db/dt = -kappa [Re(omega_211) 2 Im(omega_211) N_211 + Re(omega_322) 2 Im(omega_322) N_322]
            + kappa (2 Re(omega_211) - Re(omega_322)) R_BH
    dj/dt = -kappa [2 Im(omega_211) N_211 + 2 * 2 Im(omega_322) N_322]

(the quantum absorbed in the first process has azimuthal number 2 * 1 - 2 = 0, so it changes the hole's mass only).
Optional GW annihilation (AXION_GW=1, AXION_GW_SCALE) as in scripts/hydrogenic_spindown.jl.

Numerics. The Leaver frequencies are tabulated once on a grid in (b, a) and interpolated with a bicubic spline (the
table is saved as results/leaver_frequency_table_alpha0.3.npz and reused; a direct Leaver solve costs 0.25 s and the
model needs many thousands of evaluations). The equations are stiff once the cloud holds the hole at the superradiant
threshold, so they are integrated with LSODA, in ln N_i.

Run from the package root (numpy, scipy, mpmath; the Leaver module is imported from $AXION_MAIN_SRC, the src/
directory of the `main` branch):
    python scripts/two_level_spindown.py <f_a/M_Pl> [initial spin, default 0.9] [final t/M0, default 1e10]
Potential: AXION_POTENTIAL=cosine (default), chiral, or free (quartic_ratio = 0: reproduces leaver_spindown.py).
Writes results/two_level_<potential>_spindown_a<spin>_alpha0.3_fa<f_a/M_Pl>[_t<final t>][_gw]_diagnostics.csv with the
columns of the leaver_spindown.py output.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline

if not 2 <= len(sys.argv) <= 4:
    sys.exit("usage: two_level_spindown.py <f_a/M_Pl> [initial spin] [final t/M0]")
decay_constant_label = sys.argv[1]
spin_label = sys.argv[2] if len(sys.argv) >= 3 else "0.9"
final_time_label = sys.argv[3] if len(sys.argv) >= 4 else None
final_time = float(final_time_label) if final_time_label else 1e10
coupling_squared = float(decay_constant_label) ** 2   # kappa = (f_a/M_Pl)^2
initial_spin = float(spin_label)
initial_mass_coupling = 0.3
seed_occupations = [3e-6, 3e-6]
modes = ({"n": 2, "l": 1, "m": 1}, {"n": 3, "l": 2, "m": 2})   # |211>, |322>
results_directory = os.path.join(os.path.dirname(__file__), "..", "results")

quark_mass_ratio = 0.48   # default of ChiralPotential in src/potential.jl
potential_name = os.environ.get("AXION_POTENTIAL", "cosine")
quartic_ratio = {"cosine": 1.0, "chiral": 1 - 3 * quark_mass_ratio / (1 + quark_mass_ratio) ** 2, "free": 0.0}[potential_name]
absorption_constant = 4e-7   # 211 x 211 -> 322 x BH (Table I)
emission_constant = 1e-8     # 322 x 322 -> 211 x infinity (Table I)

gw_emission = os.environ.get("AXION_GW", "0") == "1"
gw_scale_label = os.environ.get("AXION_GW_SCALE", "1")
gw_scale = float(gw_scale_label) if gw_emission else 0.0
gw_power_table = np.loadtxt(os.path.join(os.path.dirname(__file__), "..", "data", "superrad_gw_power.csv"),
                            delimiter=",", skiprows=1)


def gw_power_coefficient(level_index, mass_coupling):
    """SuperRad coefficient of level 0 (|211>) or 1 (|322>): ln(coefficient) linear in ln alpha between table rows."""
    return np.exp(np.interp(np.log(mass_coupling), np.log(gw_power_table[:, 0]), np.log(gw_power_table[:, 1 + level_index])))


# --- Leaver frequencies M omega(alpha0 b, a) of both modes on a grid, built once ---------------------------------------
table_path = os.path.join(results_directory, "leaver_frequency_table_alpha0.3.npz")
if not os.path.exists(table_path):
    sys.path.insert(0, os.environ.get("AXION_MAIN_SRC", "/global/u1/c/clwelch/axion-sr-sim/src"))
    from common import leaver  # noqa: E402
    mass_grid = np.linspace(0.90, 1.0, 21)
    spin_grid = np.concatenate([np.linspace(0.78, 0.98, 41), np.linspace(0.9825, 0.995, 6)])
    start = time.time()
    frequencies = np.array([[[complex(leaver.solve(initial_mass_coupling * mass, spin, l=mode["l"], m=mode["m"], n=mode["n"],
                                                   spheroidal=True)) for spin in spin_grid] for mass in mass_grid]
                            for mode in modes])
    np.savez(table_path, mass_grid=mass_grid, spin_grid=spin_grid, frequencies=frequencies)
    print(f"wrote {table_path}: {frequencies.size} Leaver solves in {time.time() - start:.0f} s")
table = np.load(table_path)
real_splines = [RectBivariateSpline(table["mass_grid"], table["spin_grid"], table["frequencies"][index].real) for index in (0, 1)]
imaginary_splines = [RectBivariateSpline(table["mass_grid"], table["spin_grid"], table["frequencies"][index].imag) for index in (0, 1)]


def mode_frequencies(mass, spin):
    """Interpolated Leaver M omega of each mode in units of the current mass (complex)."""
    inside = table["mass_grid"][0] <= mass <= table["mass_grid"][-1] and table["spin_grid"][0] <= spin <= table["spin_grid"][-1]
    if not inside:
        raise ValueError(f"(b, a) = ({mass}, {spin}) is outside the Leaver frequency table")
    return [complex(real_splines[index](mass, spin)[0, 0], imaginary_splines[index](mass, spin)[0, 0]) for index in (0, 1)]


if os.environ.get("AXION_CHECK_TABLE", "0") == "1":   # spline against direct Leaver solves at points off the grid
    sys.path.insert(0, os.environ.get("AXION_MAIN_SRC", "/global/u1/c/clwelch/axion-sr-sim/src"))
    from common import leaver  # noqa: E402
    random = np.random.default_rng(1)
    worst = [0.0, 0.0, 0.0, 0.0]
    for mass, spin in zip(random.uniform(0.905, 0.999, 12), random.uniform(0.80, 0.99, 12)):
        for index, mode in enumerate(modes):
            exact = complex(leaver.solve(initial_mass_coupling * mass, spin, l=mode["l"], m=mode["m"], n=mode["n"], spheroidal=True))
            interpolated = mode_frequencies(mass, spin)[index]
            worst[index] = max(worst[index], abs(interpolated.real / exact.real - 1))
            # the rate crosses zero at the threshold, so its error is compared with the largest rate in the table
            worst[2 + index] = max(worst[2 + index], abs(interpolated.imag - exact.imag) / np.abs(table["frequencies"][index].imag).max())
    print("table check at 12 random points: largest relative error of Re omega (211, 322): %.1e, %.1e; largest error of "
          "Im omega relative to the largest rate in the table: %.1e, %.1e" % tuple(worst))


def rates(_, state):
    occupations, mass, angular_momentum = np.exp(state[:2]), state[2], state[3]
    spin, mass_coupling = angular_momentum / mass**2, initial_mass_coupling * mass
    frequencies = [frequency / mass for frequency in mode_frequencies(mass, spin)]   # M0 units
    superradiant_rates = [2 * frequency.imag * occupation for frequency, occupation in zip(frequencies, occupations)]
    absorption_rate = (quartic_ratio**2 * absorption_constant * mass_coupling**11 * (1 + np.sqrt(1 - spin**2))
                       * initial_mass_coupling * occupations[0] ** 2 * occupations[1] / mass**4)
    emission_rate = (quartic_ratio**2 * emission_constant * mass_coupling**8 * initial_mass_coupling
                     * occupations[1] ** 2 * occupations[0] / mass**4)
    gw_powers = [0.0, 0.0]
    if gw_scale != 0:
        gw_powers = [gw_scale * gw_power_coefficient(index, mass_coupling) * (coupling_squared * frequency.real * occupation / mass) ** 2
                     for index, (frequency, occupation) in enumerate(zip(frequencies, occupations))]
    gw_rates = [power / (coupling_squared * frequency.real) for power, frequency in zip(gw_powers, frequencies)]
    occupation_rates = [superradiant_rates[0] - 2 * absorption_rate + emission_rate - gw_rates[0],
                        superradiant_rates[1] + absorption_rate - 2 * emission_rate - gw_rates[1]]
    mass_rate = coupling_squared * (-sum(frequency.real * rate for frequency, rate in zip(frequencies, superradiant_rates))
                                    + (2 * frequencies[0].real - frequencies[1].real) * absorption_rate)
    angular_momentum_rate = -coupling_squared * sum(mode["m"] * rate for mode, rate in zip(modes, superradiant_rates))
    return [occupation_rates[0] / occupations[0], occupation_rates[1] / occupations[1], mass_rate, angular_momentum_rate]


saved_times = np.unique(np.concatenate([np.linspace(0, min(final_time, 2e9), 801),
                                        np.geomspace(2e9, final_time, 400)[1:] if final_time > 2e9 else []]))
start = time.time()
solution = solve_ivp(rates, (0.0, final_time), [*np.log(seed_occupations), 1.0, initial_spin], method="LSODA",
                     t_eval=saved_times, rtol=1e-9, atol=1e-11)
print(f"{potential_name} (quartic ratio {quartic_ratio:.4f}): {solution.message} {solution.nfev} rhs evaluations in "
      f"{time.time() - start:.1f} s, {len(solution.t)} saved states")

output_path = os.path.join(results_directory, f"two_level_{potential_name}_spindown_a{spin_label}_alpha0.3_fa{decay_constant_label}"
                           + (f"_t{final_time_label}" if final_time_label else "")
                           + (("_gw" + ("" if gw_scale_label == "1" else gw_scale_label)) if gw_emission else "") + "_diagnostics.csv")
with open(output_path, "w") as output:
    output.write("time,N1,N2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,b,a\n")
    for index, saved_time in enumerate(solution.t):
        log_occupation_211, log_occupation_322, mass, angular_momentum = solution.y[:, index]
        frequencies = mode_frequencies(mass, angular_momentum / mass**2)   # current-mass units, as in the other diagnostics
        row = (saved_time, np.exp(log_occupation_211), np.exp(log_occupation_322), frequencies[0].real, frequencies[1].real,
               frequencies[0].imag, frequencies[1].imag, mass, angular_momentum / mass**2)
        output.write(",".join(f"{value:.10e}" for value in row) + "\n")
final = solution.y[:, -1]
peak_index = np.argmax(solution.y[0])
print(f"peak N211 = {np.exp(solution.y[0, peak_index]):.5g} at t = {solution.t[peak_index]:.4e}; final: t = {solution.t[-1]:.3e}, "
      f"N = ({np.exp(final[0]):.6g}, {np.exp(final[1]):.6g}), a = {final[3] / final[2]**2:.5f}, M/M0 = {final[2]:.6f}\nwrote {output_path}")
