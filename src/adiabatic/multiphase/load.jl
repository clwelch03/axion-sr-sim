# Load multiphase once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_MULTIPHASE_LOADED)
    include(joinpath(@__DIR__, "../core/load.jl"))
    include(joinpath(@__DIR__, "components.jl"))
    include(joinpath(@__DIR__, "profiles.jl"))
    include(joinpath(@__DIR__, "projection.jl"))
    const ADIABATIC_MULTIPHASE_LOADED = true
end
