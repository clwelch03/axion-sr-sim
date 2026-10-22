# Load multiphase once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_MULTIPHASE_LOADED)
    include(joinpath(@__DIR__, "../core/load.jl"))
    include(joinpath(@__DIR__, "components.jl"))
    include(joinpath(@__DIR__, "radial.jl"))
    include(joinpath(@__DIR__, "profiles.jl"))
    include(joinpath(@__DIR__, "projection.jl"))
    include(joinpath(@__DIR__, "observables.jl"))
    include(joinpath(@__DIR__, "response.jl"))
    include(joinpath(@__DIR__, "evolution.jl"))
    const ADIABATIC_MULTIPHASE_LOADED = true
end
