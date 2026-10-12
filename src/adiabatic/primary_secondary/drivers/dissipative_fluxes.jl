# julia --project=. src/adiabatic/primary_secondary/drivers/dissipative_fluxes.jl primary_branch.jls secondary_branch.jls secondary_energy.csv output.csv
# The m = 0 field at ω0 = 2 Re ω1 - Re ω2 (Omiya et al. 2211.01949 Eqs. (51)-(52);
# docs/multimode-piece4-guide.md) at every saved branch point: horizon and outer
# fluxes per unit (2,2) outer amplitude of the secondary squared, the horizon
# flux per unit secondary energy F_H / Ê2, the power identity, the column-scaled
# condition number of the matching matrix, and each l's share of F_H.
# Set AXION_MAX_POINTS to stop after that many points.
using Serialization, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 4 ||
        error("require primary_branch.jls secondary_branch.jls secondary_energy.csv output.csv")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    dissipative_fluxes(ARGS...; maximum_points)
end
