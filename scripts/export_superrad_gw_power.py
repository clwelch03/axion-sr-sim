"""Export the gravitational-wave annihilation power of a scalar cloud from SuperRad to a CSV table.

SuperRad (Siemonsen, May & East, arXiv 2211.03845; Python package `superrad`) gives the power radiated by a cloud in
the level with azimuthal number m as

    P = power_coefficient_m(α) (M_cloud / M)²        (G = c = 1, so P is dimensionless),

with α = μM. Azimuthal number 1 is the |211⟩ level and azimuthal number 2 is the |322⟩ level. The coefficient is
SuperRad's `RelScalar.power_gw(m, α, spin)`:
  - m = 1: for α ≥ 0.2 a cubic interpolation of 40 Teukolsky values on 0.1756 ≤ α ≤ 0.45; for α < 0.2 the fit
    0.0109289 α¹⁴ − 0.0290106 α¹⁵;
  - m = 2: for α ≥ 0.34 a cubic interpolation of 40 Teukolsky values on 0.34 ≤ α ≤ 0.98; for α < 0.34 the fit
    6.46575e-7 α¹⁸ − 1.12205e-6 α¹⁹.
The spin argument is ignored: the values are for the spin at which the level sits at the superradiant threshold.
SuperRad's two branches do not join exactly: the coefficient jumps by +6.5% at α = 0.2 (m = 1) and by −5.2% at
α = 0.34 (m = 2). The table reproduces SuperRad as it is, including these jumps.

The table is on a uniform grid 0.010 ≤ α ≤ 0.449 with spacing 0.001 (the m = 1 data end at α = 0.45). The grid
is dense enough that interpolating ln(coefficient) linearly in ln α between rows reproduces SuperRad to the accuracy
printed at the end of this script (the two jumps excepted), so the Julia code needs no interpolation package.

Needs superrad 0.3.0 (pip install superrad==0.3.0), which is not part of the BHsuper environment.
Run from the package root:
    python scripts/export_superrad_gw_power.py [output, default data/superrad_gw_power.csv]
"""
import os
import sys
from importlib.metadata import version

import numpy as np
from superrad.rel_sca_cloud import RelScalar

output_path = sys.argv[1] if len(sys.argv) > 1 else "data/superrad_gw_power.csv"
branch_switch = {1: 0.2, 2: 0.34}  # α where SuperRad changes from the small-α fit to the interpolated data

cloud_model = RelScalar()
mass_couplings = np.round(np.arange(0.010, 0.4495, 0.001), 3)
unused_spin = 0.9  # power_gw ignores its spin argument
power_coefficients = {m: np.array([float(cloud_model.power_gw(m, alpha, unused_spin)) for alpha in mass_couplings])
                      for m in (1, 2)}

os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
with open(output_path, "w") as output:
    output.write("mass_coupling,power_coefficient_211,power_coefficient_322\n")
    for row in zip(mass_couplings, power_coefficients[1], power_coefficients[2]):
        output.write("%.3f,%.16e,%.16e\n" % row)
print(f"superrad {version('superrad')}: "
      f"wrote {len(mass_couplings)} rows to {output_path}")

# Checks. (1) the value quoted in docs/gw-emission-recommendation-2026-11-13.md; (2) the small-α fit;
# (3) log-log linear interpolation of the table against SuperRad at the midpoints between rows.
print("check 1: |211> coefficient at α = 0.3: %.4e (document: 8.4e-11)" % power_coefficients[1][mass_couplings == 0.3][0])
alpha = 0.1
print("check 2: |211> at α = 0.1: table %.6e, fit 0.0109289 α^14 − 0.0290106 α^15 = %.6e"
      % (power_coefficients[1][mass_couplings == alpha][0], 0.0109289473739731 * alpha**14 - 0.0290105840870259 * alpha**15))
midpoints = 0.5 * (mass_couplings[:-1] + mass_couplings[1:])
for m in (1, 2):
    interpolated = np.exp(np.interp(np.log(midpoints), np.log(mass_couplings), np.log(power_coefficients[m])))
    exact = np.array([float(cloud_model.power_gw(m, alpha, unused_spin)) for alpha in midpoints])
    relative_error = np.abs(interpolated / exact - 1)
    across_switch = np.abs(midpoints - branch_switch[m]) < 0.001
    print("check 3: m = %d: largest relative interpolation error at midpoints %.1e (%.1e in the interval containing "
          "the jump at α = %.2f)" % (m, relative_error[~across_switch].max(), relative_error[across_switch].max(),
                                    branch_switch[m]))
