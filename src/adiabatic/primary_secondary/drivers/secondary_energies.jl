# julia --project=. src/adiabatic/primary_secondary/drivers/secondary_energies.jl primary_branch.jls secondary_branch.jls output.csv
# Energy E2 (per unit (2,2) outer amplitude), boundary fluxes, and the energy
# balance -F/(2 Γ2 E2) of the averaged secondary at every saved branch point.
# Set AXION_MAX_POINTS to stop after that many points.
using Serialization, Printf
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 3 || error("require primary_branch.jls secondary_branch.jls output.csv")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    secondary_energies(ARGS...; maximum_points)
end
