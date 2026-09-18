@testset "Complex Newton solve" begin
    # Known complex root, with a second equation coupling the two unknowns.
    residual = state -> [state[1]^2 - (1 + 2im), state[2] - (3 - im) * state[1]]
    initial = ComplexF64[0.7 + 0.8im, 0]
    saved_initial = copy(initial)
    result = complex_newton(residual, initial)
    @test result.state ≈ [sqrt(1 + 2im), (3 - im) * sqrt(1 + 2im)]
    @test norm(residual(result.state)) < 1e-10
    @test initial == saved_initial

    # The full first Newton step from 0.1 overshoots badly; damping must help.
    damped = complex_newton(state -> state.^2 .- 1, [0.1])
    @test only(damped.state) ≈ 1
    @test_throws ErrorException complex_newton(residual, initial; max_iterations=0)
    @test_throws ArgumentError complex_newton(residual, initial; admissible=state -> false)
end

@testset "Real coordinates for complex unknowns" begin
    state = ComplexF64[1 + 2im, -3 + 4im]
    @test pack_complex(state) == [1, -3, 2, 4]
    @test unpack_complex(pack_complex(state)) == state
    @test norm(pack_complex(state)) ≈ norm(state)
    @test_throws DimensionMismatch unpack_complex([1.0, 2.0, 3.0])

    # The generic engine must keep real inputs real.
    real_result = newton_solve(values -> values.^2 .- 2, [1.0])
    @test eltype(real_result.state) <: Real
    @test only(real_result.state) ≈ sqrt(2)
end

@testset "Real Jacobian for nonanalytic complex residuals" begin
    # z + conj(z)/2 = 3+2i means 1.5x=3 and 0.5y=2.
    conjugate_residual = state -> [state[1] + 0.5conj(state[1]) - (3 + 2im)]
    result = real_newton(conjugate_residual, [0.0im])
    @test result.state ≈ [2 + 4im]
    @test norm(conjugate_residual(result.state)) < 1e-10

    # Coupling two unknowns checks the ordering of both packed halves.
    coupled_residual = state -> [
        state[1] + 0.5conj(state[1]) + state[2] - (2 + 2.5im),
        state[2] - 0.25conj(state[1]) - (-1.5 + 1.5im),
    ]
    coupled = real_newton(coupled_residual, zeros(ComplexF64, 2))
    @test coupled.state ≈ [2 + 4im, -1 + 0.5im]
    @test norm(coupled_residual(coupled.state)) < 1e-10

    # A cubic |z|^2*z term has the same conjugate dependence as our source.
    # Its unique root here is z=2+i. The callback must see complex unknowns.
    cubic_residual = state -> [state[1] * (1 + 0.2abs2(state[1])) - (4 + 2im)]
    admissible = state -> eltype(state) <: Complex && imag(state[1]) > 0
    initial = [0.1 + 0.2im]
    saved_initial = copy(initial)
    cubic = real_newton(cubic_residual, initial; admissible)
    @test cubic.state ≈ [2 + 1im]
    @test norm(cubic_residual(cubic.state)) < 1e-10
    @test cubic.residual ≈ cubic_residual(cubic.state)
    @test initial == saved_initial
    @test_throws ArgumentError real_newton(cubic_residual, [0.0im]; admissible)
    @test_throws ErrorException real_newton(cubic_residual, initial; max_iterations=0)
end

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
