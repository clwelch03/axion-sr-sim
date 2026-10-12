# Load single_cloud once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_SINGLE_CLOUD_LOADED)
    include(joinpath(@__DIR__, "../core/load.jl"))
    include(joinpath(@__DIR__, "profiles.jl"))
    include(joinpath(@__DIR__, "projection.jl"))
    include(joinpath(@__DIR__, "observables.jl"))
    include(joinpath(@__DIR__, "interpolation.jl"))
    include(joinpath(@__DIR__, "evolution.jl"))
    include(joinpath(@__DIR__, "branch_fitting.jl"))
    include(joinpath(@__DIR__, "continuation_predictor.jl"))
    include(joinpath(@__DIR__, "continuation.jl"))
    include(joinpath(@__DIR__, "branch_observables.jl"))
    include(joinpath(@__DIR__, "linear_benchmarks.jl"))
    const ADIABATIC_SINGLE_CLOUD_LOADED = true
end
