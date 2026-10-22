# The multiphase runner loads solve_two_phase_mode from the library.
# Load the single-cloud solver only to construct dilute starting guesses.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "Two-phase Newton matching prototype" begin
    grid = angular_grid(8, 9)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar, grid)
    cloud_m, spin, alpha = [1, 2], 0.99, 0.42
    components = [fourier_component(q, [l], cloud_m) for (l, q) in
        [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]]
    options = (; inner=-30.0, outer=80.0, match=5.0, step=0.4)
    # Deliberately modest resolution: test the matching algorithm, not physical
    # endpoint/truncation convergence or accurate growth rates.
    first_mode = solve_linear_mode(0.4088+1e-7im, spin, alpha, [1], 1, couplings;
        options..., tolerance=1e-9)
    second_mode = solve_linear_mode(0.4159+1e-9im, spin, alpha, [2], 2, couplings;
        options..., tolerance=1e-9)
    linear_frequencies = [first_mode.frequency, second_mode.frequency]
    frequency_guess = linear_frequencies .+ [2e-6, -2e-6]
    normalizations = [1e-4, 1e-5]
    horizon = ComplexF64[normalizations[1]*first_mode.horizon_amplitudes[1],
        normalizations[2]*second_mode.horizon_amplitudes[1], 0, 0]
    # These two outer entries must be replaced by the prescribed normalizations.
    outer = ComplexF64[99, -99, 0, 0]
    original = (copy(frequency_guess), copy(horizon), copy(outer))
    solve_trial = (frequencies, norms, cs, hs, bs; kwargs...) -> solve_two_phase_mode(
        frequencies, norms, spin, alpha, cs, cloud_m, couplings, polar, grid,
        grid.azimuth_nodes, grid.azimuth_nodes; horizon_guess=hs, outer_guess=bs,
        potential=:quartic, options..., kwargs...)
    solution = solve_trial(frequency_guess, normalizations, components, horizon, outer;
        tolerance=1e-8, max_iterations=6,
        progress_callback=(iteration, residual) -> begin
            println("Joint Newton iteration ", iteration, ": scaled residual=", residual)
            flush(stdout)
        end)
    @test solution.outer_coefficients[1:2] == normalizations
    @test solution.iterations > 0
    @test solution.scaled_residual_norm < 1e-8
    @test frequency_guess == original[1] && horizon == original[2] && outer == original[3]
    (; free_outer_coefficient_indices) = two_phase_coefficient_indices(components)
    problem = TwoPhaseProblem(; spin, alpha, components, cloud_m, couplings,
        polar_values=polar, grid, phase1_nodes=grid.azimuth_nodes,
        phase2_nodes=grid.azimuth_nodes, potential=:quartic,
        coefficient_scale=minimum(normalizations), free_outer_coefficient_indices,
        initial_bound=BitVector([real(sum(component.q .* frequency_guess)) < alpha
            for component in components]),
        integration=(; horizon_order=3, outer_order=7, options...))
    profile = TwoPhaseProfile(solution.frequencies,
        solution.horizon_coefficients, solution.outer_coefficients)
    checked = two_phase_matching_residual(profile, problem)
    @test norm(checked)/minimum(normalizations) < 1e-8
    @test isapprox(solution.residual, checked; atol=1e-18, rtol=1e-8)
    @test isapprox(solution.residual_norm, minimum(normalizations)*solution.scaled_residual_norm)
    @test maximum(abs.(solution.frequencies - linear_frequencies)) < 1e-5
    @test all(abs.(solution.outer_coefficients[3:4]) .> 1e-16)

    # Normalization is attached to the fundamental q, not its list position.
    permutation = [3, 2, 4, 1]
    reordered = solve_trial(solution.frequencies, normalizations, components[permutation],
        solution.horizon_coefficients[permutation], solution.outer_coefficients[permutation];
        tolerance=1e-8, max_iterations=0)
    @test reordered.outer_coefficients[[4, 2]] == normalizations
    @test reordered.iterations == 0
    @test reordered.horizon_coefficients ≈ solution.horizon_coefficients[permutation]
    @test reordered.outer_coefficients ≈ solution.outer_coefficients[permutation]

    @test_throws DimensionMismatch solve_trial(frequency_guess[1:1], normalizations,
        components, horizon, outer)
    @test_throws ArgumentError solve_trial(frequency_guess, [1e-4, 0.0], components, horizon, outer)
    @test_throws ArgumentError solve_trial(frequency_guess, [1e-4, 1e-5im], components, horizon, outer)
    @test_throws DimensionMismatch solve_trial(frequency_guess, normalizations,
        components, horizon[1:3], outer)
    @test_throws ArgumentError solve_trial(frequency_guess, normalizations,
        components[[1, 3, 4]], horizon[[1, 3, 4]], outer[[1, 3, 4]])
    @test_throws ArgumentError solve_trial([0.43, 0.415], normalizations, components, horizon, outer)
    @test_throws ArgumentError solve_trial([0.1, 0.3], normalizations, components, horizon, outer)
    println("Matched prototype frequencies: ", solution.frequencies)
    println("Driven outer coefficients: ", solution.outer_coefficients[3:4])
end
