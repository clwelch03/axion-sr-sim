using Test

include(joinpath(@__DIR__, "../../src/adiabatic/angular.jl"))

@testset "Angular quadrature and polar orthonormality" begin
    grid = angular_grid(16, 31)
    @test sum(grid.polar_weights) ≈ 2
    @test length(grid.azimuth_nodes) * grid.azimuth_weight ≈ 1

    polar_values = build_polar_values(5, grid.cosine_theta_nodes)
    for m in 0:5, l in m:5, other_l in m:5
        overlap = sum(
            grid.polar_weights .* polar_values[(l, m)] .* polar_values[(other_l, m)]
        )
        expected = l == other_l ? 1.0 : 0.0
        @test isapprox(overlap, expected; atol=1e-12, rtol=0)
    end
end
