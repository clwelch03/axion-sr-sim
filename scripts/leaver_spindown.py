"""Free-field spin-down with exact linear (Leaver) frequencies: the linear benchmark for the free time evolution.

Same ODE and units as scripts/hydrogenic_spindown.jl, with the hydrogenic rate and frequency replaced by the
spheroidal Leaver continued fraction of `main` (src/common/leaver.py, `solve(..., spheroidal=True)`), solved at the
current (α, a) at every right-hand-side evaluation (≈ 0.25 s per mode, cold start from the hydrogenic guess):

    d𝒩_i/dt = 2 Im ω_i 𝒩_i,   dM/dt = −κ Σ_i Re ω_i d𝒩_i/dt,   dJ/dt = −κ Σ_i m_i d𝒩_i/dt,

time in M₀, 𝒩 in f_a² M₀², b = M/M₀, j = J/M₀², α = α₀ b, a = j/b², κ = (f_a/M_Pl)². ω in M₀ units is
(M ω)/b. For a single-frequency mode the horizon fluxes are F_E = ω F_𝒩 and F_J = m F_𝒩, and the number of quanta is
an adiabatic invariant of the slowly changing background, so this is the exact linear evolution in the adiabatic
limit. Our free run (AXION_POTENTIAL=free scripts/evolution_spindown.jl) should reproduce it.

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


def mode_frequencies(mass, angular_momentum):
    """Leaver ω of each mode in M₀ units at the current background."""
    mass_coupling, spin = initial_mass_coupling * mass, angular_momentum / mass**2
    return [leaver.solve(mass_coupling, spin, l=mode["l"], m=mode["m"], n=mode["n"], spheroidal=True) / mass
            for mode in modes]


def rates(_, state):
    global evaluation_count
    evaluation_count += 1
    occupations, mass, angular_momentum = state[:2], state[2], state[3]
    frequencies = mode_frequencies(mass, angular_momentum)
    occupation_rates = [2 * frequency.imag * occupation for frequency, occupation in zip(frequencies, occupations)]
    mass_rate = -coupling_squared * sum(frequency.real * rate for frequency, rate in zip(frequencies, occupation_rates))
    angular_momentum_rate = -coupling_squared * sum(mode["m"] * rate for mode, rate in zip(modes, occupation_rates))
    return [*occupation_rates, mass_rate, angular_momentum_rate]


# saved times: dense linear grid up to 2e9 (where the free |211> stage happens at α = 0.3), then logarithmic
saved_times = np.unique(np.concatenate([np.linspace(0, min(final_time, 2e9), 401),
                                        np.geomspace(2e9, final_time, 60)[1:] if final_time > 2e9 else []]))
start = time.time()
solution = solve_ivp(rates, (0.0, final_time), [*seed_occupations, 1.0, initial_spin], method="DOP853",
                     t_eval=saved_times, rtol=1e-9, atol=[1e-14, 1e-14, 1e-13, 1e-13])
print(f"{solution.message} {evaluation_count} rhs evaluations in {time.time() - start:.0f} s, {len(solution.t)} saved states")

output_path = os.path.join(os.path.dirname(__file__), "..", "results",
                           f"leaver_spindown_a{spin_label}_alpha0.3_fa{decay_constant_label}_diagnostics.csv")
with open(output_path, "w") as output:
    output.write("time,N1,N2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,b,a\n")
    for index, saved_time in enumerate(solution.t):
        occupation_211, occupation_322, mass, angular_momentum = solution.y[:, index]
        # current-mass units, as in the evolution_spindown.jl diagnostics
        frequencies = [frequency * mass for frequency in mode_frequencies(mass, angular_momentum)]
        row = (saved_time, occupation_211, occupation_322, frequencies[0].real, frequencies[1].real,
               frequencies[0].imag, frequencies[1].imag, mass, angular_momentum / mass**2)
        output.write(",".join(f"{value:.10e}" for value in row) + "\n")
final = solution.y[:, -1]
print(f"final: t = {solution.t[-1]:.3e}, N = ({final[0]:.6g}, {final[1]:.6g}), a = {final[3] / final[2]**2:.5f}, "
      f"M/M0 = {final[2]:.6f}\nwrote {output_path}")
