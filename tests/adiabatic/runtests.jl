using Test

include("test_loaders.jl")

@testset "core" begin
    include("core/runtests.jl")
end
@testset "single_cloud" begin
    include("single_cloud/runtests.jl")
end
@testset "primary_secondary" begin
    include("primary_secondary/runtests.jl")
end
