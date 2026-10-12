# Load core once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_CORE_LOADED)
    include(joinpath(@__DIR__, "angular.jl"))
    include(joinpath(@__DIR__, "geometry.jl"))
    include(joinpath(@__DIR__, "potentials.jl"))
    include(joinpath(@__DIR__, "radial_operator.jl"))
    include(joinpath(@__DIR__, "boundary_series.jl"))
    include(joinpath(@__DIR__, "boundaries.jl"))
    include(joinpath(@__DIR__, "integration.jl"))
    include(joinpath(@__DIR__, "newton.jl"))
    include(joinpath(@__DIR__, "quadrature.jl"))
    include(joinpath(@__DIR__, "sensitivity.jl"))
    include(joinpath(@__DIR__, "units.jl"))
    const ADIABATIC_CORE_LOADED = true
end
