"""Free-field spin-down with exact linear (Leaver) frequencies: the linear benchmark for the free time evolution.

Same ODE and units as scripts/hydrogenic_spindown.jl, with the hydrogenic rate and frequency replaced by the
spheroidal Leaver continued fraction of `main` (src/common/leaver.py, `solve(..., spheroidal=True)`), solved at the
current (α, a) at every right-hand-side evaluation (≈ 0.25 s per mode, cold start from the hydrogenic guess):

    d𝒩_i/dt = 2 Im ω_i 𝒩_i,   dM/dt = −κ Σ_i Re ω_i d𝒩_i/dt,   dJ/dt = −κ Σ_i m_i d𝒩_i/dt,

time in M₀, 𝒩 in f_a² M₀², b = M/M₀, j = J/M₀², α = α₀ b, a = j/b², κ = (f_a/M_Pl)². ω in M₀ units is
(M ω)/b. For a single-frequency mode the horizon fluxes are F_E = ω F_𝒩 and F_J = m F_𝒩, and the number of quanta is
an adiabatic invariant of the slowly changing background, so this is the exact linear evolution in the adiabatic
limit. Our free run (AXION_POTENTIAL=free scripts/evolution_spindown.jl) should reproduce it.

Optional gravitational-wave (GW) annihilation, AXION_GW=1 (default off), exactly as in scripts/hydrogenic_spindown.jl
(see the header there): each level also loses quanta at the rate scale × power_coefficient_i(α) κ Re ω_i 𝒩_i² / b²,
with the SuperRad coefficient from data/superrad_gw_power.csv and scale = AXION_GW_SCALE (default 1). The hole
equations keep only the superradiant term. The radiated energy E_GW/M₀ is a fifth state variable. The output name
then ends in _gw_diagnostics.csv (_gw<scale>_diagnostics.csv for a scale other than 1).

Optional perturbative self-interaction, AXION_QUARTIC_COUPLING=c (default 0: none). This turns the script into the
two-level model of Baryakhtar, Galanis, Lasenby & Simon, arXiv 2011.11646 (docs/papers/2011.11646.pdf), Eqs. (34)-(37),
keeping the two processes they single out, with the leading-order rate constants of their Table I:
  - two |211> quanta -> one |322> quantum + one quantum absorbed by the hole, rate constant 4e-7 α¹¹ (1 + √(1 − a²)) μ;
  - two |322> quanta -> one |211> quantum + one quantum radiated to infinity, rate constant 1e-8 α⁸ μ.
Both are second order in the quartic coupling λ of the potential. c is λ in units of the cosine potential's
λ = μ²/f_a² (c = 1: cosine; c = 0.3426: the quartic term of the chiral potential with quark mass ratio 0.48), so the
rate constants are multiplied by c². In this script's occupations 𝒩 = N/(f_a² M₀²) (N the number of quanta), and with
the paper's normalized occupation ε = N/(G M²) = κ 𝒩 / b², the terms are
    d𝒩₂₁₁/dt = … − 2 R_BH 𝒩₂₁₁² 𝒩₃₂₂ + R_∞ 𝒩₃₂₂² 𝒩₂₁₁,      d𝒩₃₂₂/dt = … + R_BH 𝒩₂₁₁² 𝒩₃₂₂ − 2 R_∞ 𝒩₃₂₂² 𝒩₂₁₁,
    R_BH = c² 4e-7 α¹¹ (1 + √(1 − a²)) α₀ / b⁴,   R_∞ = c² 1e-8 α⁸ α₀ / b⁴     (per unit t/M₀; R_BH and R_∞ are names
used only here). f_a drops out of these terms, as in the full code. The quantum absorbed by the hole in the first
process has energy 2 Re ω₂₁₁ − Re ω₃₂₂ and azimuthal number 2·1 − 2 = 0, so it adds κ (2 Re ω₂₁₁ − Re ω₃₂₂) R_BH 𝒩₂₁₁² 𝒩₃₂₂
to db/dt and nothing to dj/dt (the paper's Eq. 37 has the same term with the energy approximated by μ). The linear
rates and frequencies remain the Leaver ones (the paper's leading-order superradiance rates are not used). The
Table I constants are leading order in α and the paper itself says they can be poor at the α of interest (its p. 11).
The output name then gets _quartic<c> before _diagnostics.csv.

Integrator: AXION_ODE_METHOD, a scipy solve_ivp method name (default DOP853, explicit). Once the |211> cloud sits at
its threshold the equations are stiff (the hole's spin relaxes to the threshold much faster than anything else
changes), and an explicit method then takes very many small steps; LSODA switches to an implicit method there.
With a method other than DOP853 the output name gets the suffix _<method> before _diagnostics.csv.

Run from the package root (needs mpmath, numpy, scipy; on Perlmutter `module load python`):
    python scripts/leaver_spindown.py <f_a/M_Pl> [initial spin, default 0.9] [final t/M₀, default 1e10]
The Leaver module is imported from $AXION_MAIN_SRC (default: the main checkout's src/).
Writes results/leaver_spindown_a<spin>_alpha0.3_fa<f_a/M_Pl>_diagnostics.csv (one row per saved time).
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.environ.get("AXION_MAIN_SRC", "/global/u1/c/clwelch/axion-sr-sim/src"))
from common import leaver  # noqa: E402

if not 2 <= len(sys.argv) <= 4:
    sys.exit("usage: leaver_spindown.py <f_a/M_Pl> [initial spin] [final t/M0]")
decay_constant_label = sys.argv[1]
spin_label = sys.argv[2] if len(sys.argv) >= 3 else "0.9"
final_time = float(sys.argv[3]) if len(sys.argv) >= 4 else 1e10
coupling_squared = float(decay_constant_label) ** 2
initial_spin = float(spin_label)
initial_mass_coupling = 0.3
seed_occupations = [3e-6, 3e-6]
modes = ({"n": 2, "l": 1, "m": 1}, {"n": 3, "l": 2, "m": 2})   # |211>, |322>
evaluation_count = 0
ode_method = os.environ.get("AXION_ODE_METHOD", "DOP853")
gw_emission = os.environ.get("AXION_GW", "0") == "1"
gw_scale_label = os.environ.get("AXION_GW_SCALE", "1")
gw_scale = float(gw_scale_label) if gw_emission else 0.0
quartic_coupling_label = os.environ.get("AXION_QUARTIC_COUPLING", "0")
quartic_coupling = float(quartic_coupling_label)
gw_power_table = np.loadtxt(os.path.join(os.path.dirname(__file__), "..", "data", "superrad_gw_power.csv"),
                            delimiter=",", skiprows=1)


def gw_power_coefficient(level_index, mass_coupling):
    """SuperRad coefficient of level 0 (|211>) or 1 (|322>): ln(coefficient) linear in ln α between table rows."""
    if not gw_power_table[0, 0] <= mass_coupling <= gw_power_table[-1, 0]:
        raise ValueError(f"α = {mass_coupling} is outside the GW power table")
    return np.exp(np.interp(np.log(mass_coupling), np.log(gw_power_table[:, 0]),
                            np.log(gw_power_table[:, 1 + level_index])))


def mode_frequencies(mass, angular_momentum):
    """Leaver ω of each mode in M₀ units at the current background."""
    mass_coupling, spin = initial_mass_coupling * mass, angular_momentum / mass**2
    return [leaver.solve(mass_coupling, spin, l=mode["l"], m=mode["m"], n=mode["n"], spheroidal=True) / mass
            for mode in modes]


def rates(_, state):
    global evaluation_count
    evaluation_count += 1
    if evaluation_count % 200 == 0:   # progress line: the script otherwise prints nothing until the end
        print(f"evaluation {evaluation_count}: t = {_:.6e}, N = ({state[0]:.6g}, {state[1]:.6g}), "
              f"{time.time() - start:.0f} s", flush=True)
    occupations, mass, angular_momentum = state[:2], state[2], state[3]
    frequencies = mode_frequencies(mass, angular_momentum)
    superradiant_rates = [2 * frequency.imag * occupation for frequency, occupation in zip(frequencies, occupations)]
    mass_rate = -coupling_squared * sum(frequency.real * rate for frequency, rate in zip(frequencies, superradiant_rates))
    angular_momentum_rate = -coupling_squared * sum(mode["m"] * rate for mode, rate in zip(modes, superradiant_rates))
    # GW annihilation: P_i = scale × coefficient_i(α) (M_cloud,i / M)², with M_cloud,i / M₀ = κ Re ω_i 𝒩_i
    gw_powers = [0.0, 0.0]
    if gw_scale != 0:
        gw_powers = [gw_scale * gw_power_coefficient(index, initial_mass_coupling * mass)
                     * (coupling_squared * frequency.real * occupation / mass) ** 2
                     for index, (frequency, occupation) in enumerate(zip(frequencies, occupations))]
    occupation_rates = [rate - power / (coupling_squared * frequency.real)
                        for rate, power, frequency in zip(superradiant_rates, gw_powers, frequencies)]
    if quartic_coupling != 0:   # Baryakhtar et al. Eqs. (34)-(37), Table I rate constants (see the header)
        mass_coupling, spin = initial_mass_coupling * mass, angular_momentum / mass**2
        common_factor = quartic_coupling**2 * initial_mass_coupling / mass**4
        to_322_and_hole = common_factor * 4e-7 * mass_coupling**11 * (1 + np.sqrt(1 - spin**2)) * occupations[0]**2 * occupations[1]
        to_211_and_infinity = common_factor * 1e-8 * mass_coupling**8 * occupations[1]**2 * occupations[0]
        occupation_rates[0] += -2 * to_322_and_hole + to_211_and_infinity
        occupation_rates[1] += to_322_and_hole - 2 * to_211_and_infinity
        mass_rate += coupling_squared * (2 * frequencies[0].real - frequencies[1].real) * to_322_and_hole
    return [*occupation_rates, mass_rate, angular_momentum_rate, sum(gw_powers)]


# saved times: dense linear grid up to 2e9 (where the free |211> stage happens at α = 0.3), then logarithmic
saved_times = np.unique(np.concatenate([np.linspace(0, min(final_time, 2e9), 401),
                                        np.geomspace(2e9, final_time, 60)[1:] if final_time > 2e9 else []]))
start = time.time()
solution = solve_ivp(rates, (0.0, final_time), [*seed_occupations, 1.0, initial_spin, 0.0], method=ode_method,
                     t_eval=saved_times, rtol=1e-9, atol=[1e-14, 1e-14, 1e-13, 1e-13, 1e-13])
print(f"{solution.message} {evaluation_count} rhs evaluations in {time.time() - start:.0f} s, {len(solution.t)} saved states")

output_path = os.path.join(os.path.dirname(__file__), "..", "results",
                           f"leaver_spindown_a{spin_label}_alpha0.3_fa{decay_constant_label}"
                           + (("_gw" + ("" if gw_scale_label == "1" else gw_scale_label)) if gw_emission else "")
                           + ("" if quartic_coupling == 0 else "_quartic" + quartic_coupling_label)
                           + ("" if ode_method == "DOP853" else "_" + ode_method.lower()) + "_diagnostics.csv")
with open(output_path, "w") as output:
    output.write("time,N1,N2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,b,a" + (",radiated_energy" if gw_emission else "") + "\n")
    for index, saved_time in enumerate(solution.t):
        occupation_211, occupation_322, mass, angular_momentum, radiated_energy = solution.y[:, index]
        # current-mass units, as in the evolution_spindown.jl diagnostics
        frequencies = [frequency * mass for frequency in mode_frequencies(mass, angular_momentum)]
        row = (saved_time, occupation_211, occupation_322, frequencies[0].real, frequencies[1].real,
               frequencies[0].imag, frequencies[1].imag, mass, angular_momentum / mass**2,
               *((radiated_energy,) if gw_emission else ()))
        output.write(",".join(f"{value:.10e}" for value in row) + "\n")
final = solution.y[:, -1]
print(f"final: t = {solution.t[-1]:.3e}, N = ({final[0]:.6g}, {final[1]:.6g}), a = {final[3] / final[2]**2:.5f}, "
      f"M/M0 = {final[2]:.6f}\nwrote {output_path}"
      + (f"\nradiated in GW: E_GW/M0 = {final[4]:.6e} (GW scale {gw_scale_label})" if gw_emission else ""))
