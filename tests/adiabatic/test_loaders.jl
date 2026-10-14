using Test

# Isolated calling modules expose accidental dependencies hidden by Main's
# already-loaded functions. Loading a model twice must not redefine methods.
# Explicit target modules keep the editor from assigning shared source files
# to the first test module that includes them, hiding them from other models.
module CoreLoaderCheck
Base.include(@__MODULE__, joinpath(@__DIR__, "../../src/adiabatic/core/load.jl"))
end
module SingleLoaderCheck
Base.include(@__MODULE__, joinpath(@__DIR__, "../../src/adiabatic/single_cloud/load.jl"))
end
module MultiphaseLoaderCheck
Base.include(@__MODULE__, joinpath(@__DIR__, "../../src/adiabatic/multiphase/load.jl"))
end
module FullLoaderCheck
Base.include(@__MODULE__, joinpath(@__DIR__, "../../src/adiabatic/load.jl"))
end

@testset "Model loader boundaries and repeated loading" begin
    @test CoreLoaderCheck.potential_value(0.0) == 0.0
    @test !isdefined(CoreLoaderCheck, :solve_nonlinear_mode)
    @test isdefined(SingleLoaderCheck, :solve_nonlinear_mode)
    @test !isdefined(SingleLoaderCheck, :solve_secondary_mode)
    @test isdefined(MultiphaseLoaderCheck, :radial_operator_coefficients)
    @test !isdefined(MultiphaseLoaderCheck, :solve_secondary_mode)
    @test isdefined(FullLoaderCheck, :solve_secondary_mode)
    @test isdefined(FullLoaderCheck, :qcd_stationary_point)
    @test !isdefined(FullLoaderCheck, :ADIABATIC_MULTIPHASE_LOADED)
    method = first(methods(FullLoaderCheck.solve_nonlinear_mode))
    Base.include(FullLoaderCheck, joinpath(@__DIR__, "../../src/adiabatic/single_cloud/load.jl"))
    Base.include(FullLoaderCheck, joinpath(@__DIR__, "../../src/adiabatic/load.jl"))
    @test first(methods(FullLoaderCheck.solve_nonlinear_mode)) === method
end
