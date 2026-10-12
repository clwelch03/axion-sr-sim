# julia --project=. src/adiabatic/primary_secondary/drivers/continue_secondary.jl primary_branch.jls output_prefix
# Solve the secondary's Eq. (61) (Omiya et al. 2211.01949) at every saved primary
# point, in amplitude order. Set AXION_MAX_POINTS to stop after that many points.
using Serialization, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("require primary_branch.jls output_prefix")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    # Mμ = 0.3 (outer r* = 300) needs AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2.
    tolerance = parse(Float64, get(ENV, "AXION_SECONDARY_TOLERANCE", "1e-8"))
    ulp_tolerance_factor = haskey(ENV, "AXION_SECONDARY_ULP_TOLERANCE_FACTOR") ?
        parse(Float64, ENV["AXION_SECONDARY_ULP_TOLERANCE_FACTOR"]) : nothing
    continue_secondary(ARGS[1], ARGS[2]; tolerance, ulp_tolerance_factor, maximum_points)
end
