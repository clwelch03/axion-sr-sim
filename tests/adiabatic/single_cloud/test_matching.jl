@testset "Bound-mode frequency and normalization" begin
    grid = angular_grid(16, 31)
    couplings = build_cosine_squared_couplings(
        5, build_polar_values(5, grid.cosine_theta_nodes), grid)
    spin, alpha, guess = 0.99, 0.42, 0.4088 + 1e-7im
    solution = solve_linear_mode(guess, spin, alpha, [1, 3, 5], 1, couplings)
    @test solution.outer_amplitudes[1] == 1
    @test 0 < real(solution.frequency) < alpha
    @test imag(solution.frequency) > 0
    @test solution.residual_norm < 1e-10
    @test 0 < solution.iterations <= 15
    checked_residual = matching_residual(solution.frequency, solution.horizon_amplitudes,
        solution.outer_amplitudes, spin, alpha, [1, 3, 5], 1, couplings)
    @test norm(checked_residual) < 1e-10

    # Smaller RK4 spacing, farther endpoints, and a different matching point.
    refined = solve_linear_mode(solution.frequency, spin, alpha, [1, 3, 5], 1, couplings;
        horizon_guess=solution.horizon_amplitudes, outer_tail_guess=solution.outer_amplitudes[2:end],
        step=0.025, inner=-140.0, outer=200.0, match=10.0, tolerance=1e-12)
    @test abs(real(refined.frequency - solution.frequency)) < 1e-11
    # Test the tiny growth rate separately: a complex relative test can hide errors in it.
    @test isapprox(imag(refined.frequency), imag(solution.frequency); rtol=1e-7)

    # Independent reference: common.leaver.solve(0.42, 0.99, spheroidal=True, N=1000).
    # That code uses Lambda=2-c^2/5, exactly our single-l=1 angular truncation.
    # Do not use this reference as an exact target for the coupled [1,3,5] solve.
    leaver = 0.40880556672489204 + 1.5038614826273717e-7im
    single = solve_linear_mode(guess, spin, alpha, [1], 1, couplings;
        step=0.025, inner=-120.0, outer=150.0, tolerance=1e-12)
    @test abs(real(single.frequency - leaver)) < 1e-11
    @test isapprox(imag(single.frequency), imag(leaver); rtol=1e-7)
end
