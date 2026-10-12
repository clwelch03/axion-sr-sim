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
