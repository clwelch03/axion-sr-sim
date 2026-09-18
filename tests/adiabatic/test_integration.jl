using LinearAlgebra: norm

@testset "RK4 radial integration" begin
    # R'=Q, Q'=-R has the exact solution [cos(x), -sin(x)].
    function oscillator!(derivative, state, position)
        derivative[1] = state[2]
        derivative[2] = -state[1]
    end
    for (start, stop) in ((0.0, 1.0), (1.0, 0.0))
        initial = [cos(start), -sin(start)]
        saved_initial = copy(initial)
        exact = [cos(stop), -sin(stop)]
        errors = [norm(integrate_radial(oscillator!, initial, start, stop; step) - exact)
                  for step in (0.1, 0.05, 0.025)]
        # Fourth-order global error falls by 2^4 when the spacing is halved.
        @test 15 < errors[1] / errors[2] < 17
        @test 15 < errors[2] / errors[3] < 17
        @test initial == saved_initial
    end

    # A position-dependent RHS checks stage positions and a shortened spacing.
    # y'=3i*x^2 integrates exactly to y=i*x^3 under RK4.
    polynomial! = (derivative, state, position) -> (derivative .= 3im * position^2)
    for (start, stop) in ((0.0, 1.0), (1.0, 0.0))
        result = integrate_radial(polynomial!, [im * start^3], start, stop; step=0.13)
        @test isapprox(only(result), im * stop^3; atol=1e-14)
    end
    @test integrate_radial(oscillator!, [1, 0], 0.0, 0.0) == [1, 0]
    @test_throws ArgumentError integrate_radial(oscillator!, [1, 0], 0.0, 1.0; step=0)
end

@testset "Linear boundary matching" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    spin, alpha, frequency = 0.99, 0.42, 0.4088 + 1.5e-7im
    l_values = [1, 3, 5]
    horizon_amplitudes = ComplexF64[1, 0.01im, -0.002]
    outer_amplitudes = ComplexF64[1, -0.02, 0.001im]
    saved_horizon, saved_outer = copy(horizon_amplitudes), copy(outer_amplitudes)
    residual = matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                                 spin, alpha, l_values, 1, couplings)
    @test length(residual) == 6
    @test all(isfinite, residual)
    @test horizon_amplitudes == saved_horizon
    @test outer_amplitudes == saved_outer

    # Reordering l must reorder BOTH residual halves without changing physics.
    permutation = [3, 1, 2]
    reordered = matching_residual(frequency, horizon_amplitudes[permutation],
        outer_amplitudes[permutation], spin, alpha, l_values[permutation], 1, couplings)
    @test isapprox(reordered, residual[vcat(permutation, permutation .+ 3)]; rtol=1e-10)
    @test_throws ArgumentError matching_residual(frequency, horizon_amplitudes,
        outer_amplitudes, spin, alpha, l_values, 1, couplings; match=100.0)
end
