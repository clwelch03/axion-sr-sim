# Gravitational-wave annihilation power of a scalar cloud, from the SuperRad table data/superrad_gw_power.csv
# (written by scripts/export_superrad_gw_power.py; see that file for where the numbers come from and their limits).
#
# A cloud of mass M_cloud in one level, around a hole of mass M, radiates the power
#     P = power_coefficient(α) (M_cloud / M)²        (G = c = 1, so P is dimensionless),
# with α = μM the current mass coupling. Level 1 is |211⟩, level 2 is |322⟩.

"""
    read_gw_power_table(path)

Read the table. Returns the logarithms of the three columns, `(log_mass_coupling, log_coefficients)`, where
`log_coefficients[level]` is the column of level 1 (|211⟩) or 2 (|322⟩).
"""
function read_gw_power_table(path=joinpath(@__DIR__, "..", "data", "superrad_gw_power.csv"))
    rows = [parse.(Float64, split(line, ',')) for line in readlines(path)[2:end]]
    (log_mass_coupling=[log(row[1]) for row in rows],
     log_coefficients=([log(row[2]) for row in rows], [log(row[3]) for row in rows]))
end

"""
    gw_power_coefficient(table, level, mass_coupling)

The coefficient of level 1 (|211⟩) or 2 (|322⟩) at `mass_coupling`. Between two rows of the table, ln(coefficient) is
taken to be a straight line in ln α. This is exact for a power law, and the coefficient is close to one (≈ α¹⁴ for
|211⟩). Outside the table (α < 0.010 or α > 0.449) it is an error.
"""
function gw_power_coefficient(table, level, mass_coupling)
    log_mass_coupling = log(mass_coupling)
    grid = table.log_mass_coupling
    grid[1] ≤ log_mass_coupling ≤ grid[end] || error("α = $mass_coupling is outside the GW power table")
    upper = clamp(searchsortedfirst(grid, log_mass_coupling), 2, length(grid))   # row at or above α
    weight = (log_mass_coupling - grid[upper - 1]) / (grid[upper] - grid[upper - 1])
    column = table.log_coefficients[level]
    exp((1 - weight) * column[upper - 1] + weight * column[upper])
end
