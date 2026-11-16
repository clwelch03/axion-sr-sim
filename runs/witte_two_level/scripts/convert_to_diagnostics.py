"""Convert the raw output of run_witte_two_level.jl to our diagnostics CSV (columns time,N1,N2,b,a,threshold).

Usage: python convert_to_diagnostics.py <raw csv> <f_a/M_Pl> <output csv> [stride]
stride (default 1) keeps every stride-th row plus the last; their output has ~1e6 rows on a log time grid
Conversions use the constants of Witte & Mummery's code (Core/constants.jl), the same ones run_witte_two_level.jl uses
for the forward conversion: M0 = G_NEW * M_BH * HBAR seconds, one of their years = 3.15e7 s. (With the physical
GM_sun/c^3 = 4.925491e-6 s and 3.15576e7 s/yr the year per M0 would be 0.2% smaller; using the code's constants in both
directions keeps t/M0 consistent with the dynamics, because the rates enter as rate_eV * G_NEW * M_BH in M0 units.)
time = t_years / years_per_M0 [M0], N = u / kappa [f_a^2 M0^2], b = M / M_BH_initial, threshold = 4 Omega / (1 + 4 Omega^2)
with Omega = Re M omega_211 (current mass, current spin) interpolated from results/leaver_frequency_table_alpha0.3.npz
by a RectBivariateSpline in (b, a), as mode_frequencies in scripts/two_level_spindown.py.
"""
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline

G_NEW, HBAR, YEAR_SECONDS, MASS_BH = 7484169213.942707, 6.582119e-16, 3.15e7, 10.0
raw_path, ratio_text, output_path = sys.argv[1:4]
stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1
kappa = float(ratio_text) ** 2
years_per_initial_mass = G_NEW * MASS_BH * HBAR / YEAR_SECONDS
raw = np.loadtxt(raw_path, delimiter=",", skiprows=1)
keep = np.unique(np.append(np.arange(0, len(raw), stride), len(raw) - 1))
raw = raw[keep]
time_years, u_211, u_322, spin, mass = raw.T
table = np.load("/global/scratch/projects/pc_heptheory/clwelch/axion-sr-sim/results/leaver_frequency_table_alpha0.3.npz")
real_spline = RectBivariateSpline(table["mass_grid"], table["spin_grid"], table["frequencies"][0].real)
b = mass / MASS_BH
omega = real_spline.ev(b, spin)
threshold = 4 * omega / (1 + 4 * omega**2)
with open(output_path, "w") as out:
    out.write("time,N1,N2,b,a,threshold\n")
    for row in zip(time_years / years_per_initial_mass, u_211 / kappa, u_322 / kappa, b, spin, threshold):
        out.write(",".join("%.10e" % value for value in row) + "\n")
print("wrote", output_path, len(time_years), "rows; years per M0 =", years_per_initial_mass)
