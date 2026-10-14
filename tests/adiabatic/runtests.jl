using Test

include("test_loaders.jl")

# Explicit targets keep macro-wrapped test execution out of the editor's
# source include tree; each group's runner provides its analysis context.
@testset "core" begin
    Base.include(@__MODULE__, joinpath(@__DIR__, "core/runtests.jl"))
end
@testset "single_cloud" begin
    Base.include(@__MODULE__, joinpath(@__DIR__, "single_cloud/runtests.jl"))
end
@testset "primary_secondary" begin
    Base.include(@__MODULE__, joinpath(@__DIR__, "primary_secondary/runtests.jl"))
end
@testset "multiphase" begin
    Base.include(@__MODULE__, joinpath(@__DIR__, "multiphase/runtests.jl"))
end
