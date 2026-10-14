using Test
include(joinpath(@__DIR__, "../../../src/adiabatic/multiphase/load.jl"))
include("test_components.jl")
include("test_radial_operator.jl")
include("test_phase_separation.jl")
include("test_nonlinear_source.jl")
include("test_source_validation.jl")
include("test_nonlinear_rhs.jl")
