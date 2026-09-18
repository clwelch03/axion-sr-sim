using Test

include(joinpath(@__DIR__, "../../src/adiabatic/angular.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/nonlinear.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/tortoise.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/radial.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/boundary_series.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/boundaries.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/integration.jl"))
include(joinpath(@__DIR__, "../../src/adiabatic/matching.jl"))

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

@testset "Angular couplings" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(5, grid.cosine_theta_nodes)

    couplings = build_cosine_squared_couplings(5, polar_values, grid)

    @test couplings[(1, 1, 1)] ≈ 1 / 5
    @test couplings[(1, 3, 1)] ≈ sqrt(8 / 175)
    @test isapprox(couplings[(1, 5, 1)], 0; atol=1e-12)
end

@testset "Radial operator coefficients" begin
    frequency = 0.4 + 0.02im

    # Schwarzschild: B = r^2 * (omega^2 / (1 - 2/r) - mu^2) - l(l+1).
    schwarzschild = radial_operator_coefficients(4.0, 0.0, 0.5, frequency, 1, 1)
    @test schwarzschild.delta == 8
    @test schwarzschild.delta_prime == 6
    @test schwarzschild.local_coefficient ≈ 32frequency^2 - 6
    @test schwarzschild.coupling_scale == 0

    # At r=2, a=1/2, mu=1/4, l=2, m=1, the Kerr operator gives
    # B = 72*omega^2 - 16*omega - 21/4, excluding the cos^2(theta) term.
    kerr = radial_operator_coefficients(2.0, 0.5, 0.25, frequency, 2, 1)
    @test kerr.delta == 1 / 4
    @test kerr.delta_prime == 2
    @test kerr.local_coefficient ≈ 72frequency^2 - 16frequency - 21 / 4
    @test kerr.coupling_scale ≈ frequency^2 / 4 - 1 / 64

    # At omega=mu, the linear angular mixing vanishes even with nonzero spin.
    threshold = radial_operator_coefficients(2.0, 0.5, 0.25, 0.25, 2, 1)
    @test threshold.coupling_scale == 0
end

include("test_tortoise.jl")
include("test_boundaries.jl")
include("test_boundary_state.jl")
include("test_integration.jl")
include("test_matching.jl")
include("test_nonlinear.jl")
include("test_nonlinear_radial.jl")
include("test_nonlinear_integration.jl")
include("test_nonlinear_matching.jl")
