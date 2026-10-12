# julia --project=. src/adiabatic/single_cloud/drivers/branch_energy_derivatives.jl \
#     checkpoint.jls output_prefix [shard_index=0] [shard_count=1]
# Recompute energy, signed fluxes, and their amplitude derivatives at saved points.
using Serialization, LinearAlgebra, SHA
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    2 <= length(ARGS) <= 4 ||
        error("require checkpoint.jls output_prefix [shard_index] [shard_count]")

    branch_energy_derivatives(
        ARGS[1], ARGS[2];
        shard_index=length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0,
        shard_count=length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1,
    )
end
