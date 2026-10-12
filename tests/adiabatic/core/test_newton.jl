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

@testset "Newton correction until convergence" begin
    # This distant initial guess needs more than ten corrections to reach sqrt(2).
    residual = state -> state.^2 .- 2
    @test_throws ErrorException newton_solve(residual, [1e5]; max_iterations=10)
    history = Tuple{Int,Float64}[]
    solution = newton_solve(residual, [1e5]; max_iterations=nothing,
        progress_callback=(iteration, mismatch) -> push!(history, (iteration, mismatch)))
    @test solution.iterations > 10
    @test only(solution.state) ≈ sqrt(2)
    @test norm(solution.residual) <= 1e-10
    @test length(history) == solution.iterations + 1
    @test all(diff(last.(history)) .< 0)

    # An intentionally strict progress requirement exercises the stagnation guard.
    @test_throws ErrorException newton_solve(residual, [1e5];
        max_iterations=nothing, minimum_progress=0.99, stagnation_limit=2)
    @test_throws ArgumentError newton_solve(residual, [1.0]; minimum_progress=0)
    @test_throws ArgumentError newton_solve(residual, [1.0]; stagnation_limit=0)
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
