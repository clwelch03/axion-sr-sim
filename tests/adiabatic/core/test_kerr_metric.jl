# Candidate checkpoint: load kerr_metric before this file.
using LinearAlgebra: I, Diagonal, det, issymmetric

@testset "Kerr covariant and inverse metrics" begin
    for spin in (-0.7, 0.0, 0.99), position in (-5.0, 5.0, 30.0), x in (-0.8, 0.0, 0.7)
        geometry = tortoise_geometry(position, spin)
        data = kerr_metric(position, spin, x)
        @test issymmetric(data.metric) && issymmetric(data.inverse_metric)
        @test isapprox(data.metric * data.inverse_metric, Matrix{Float64}(I, 4, 4);
            atol=2e-12, rtol=2e-12)
        @test isapprox(det(data.metric),
            -(geometry.radius^2 + spin^2*x^2)^2 * (1-x^2); rtol=2e-12)
    end
    position, x = 10.0, 0.3
    radius = tortoise_geometry(position, 0.0).radius
    lapse_squared = 1 - 2/radius
    expected = Diagonal([-lapse_squared, 1/lapse_squared, radius^2, radius^2*(1-x^2)])
    schwarzschild = kerr_metric(position, 0.0, x)
    @test schwarzschild.metric ≈ expected
    @test schwarzschild.inverse_metric ≈ inv(expected)
    # The cross coefficient of dt*dvarphi in ds^2 is twice g_tvarphi.
    spin = 0.7
    radius = tortoise_geometry(position, spin).radius
    @test 2kerr_metric(position, spin, x).metric[1, 4] ≈
        -4spin*radius*(1-x^2)/(radius^2 + spin^2*x^2)
    @test_throws DomainError kerr_metric(position, spin, 1.0)
    @test_throws DomainError kerr_metric(position, spin, NaN)
end
