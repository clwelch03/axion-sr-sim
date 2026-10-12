# Load multiphase once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_MULTIPHASE_LOADED)
    include(joinpath(@__DIR__, "../core/load.jl"))
    const ADIABATIC_MULTIPHASE_LOADED = true
end
