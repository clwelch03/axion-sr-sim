# julia --project=. src/adiabatic/primary_secondary/drivers/radiative_fluxes.jl primary_branch.jls secondary_branch.jls secondary_energy.csv output.csv
# The m = 3 field at ω3 = 2 Re ω2 - Re ω1 (Omiya et al. 2211.01949 Eqs. (51)-(52))
# at every saved branch point: fluxes to infinity and into the horizon per unit
# (2,2) outer amplitude of the secondary to the fourth power, the same divided by
# Ê2² (flux per unit secondary energy squared), the power identity, the
# column-scaled condition number of the matching matrix, and each l's share of
# the flux to infinity. Set AXION_MAX_POINTS to stop after that many points and
# AXION_RADIATIVE_MAX_L (odd, default 7) to change the l truncation.
using Serialization, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 4 ||
        error("require primary_branch.jls secondary_branch.jls secondary_energy.csv output.csv")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    maximum_l = parse(Int, get(ENV, "AXION_RADIATIVE_MAX_L", "7"))
    radiative_fluxes(ARGS...; maximum_l, maximum_points)
end
