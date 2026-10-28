# Small, deterministic regression for the verification machinery and both
# derivative routes through a short nonlinear matching calculation.
using Test
include(joinpath(@__DIR__, "../../../src/adiabatic/multiphase/autodiff_vs_optimized_test.jl"))

@testset "AD comparison and saved diagnostics" begin
    toy = state -> [state[1]^2 + 3state[2], state[1]*state[2]]
    state = [2.0, 5.0]
    mktempdir() do directory
        path = joinpath(directory, "comparison.jls")
        comparison = compare_ad_jacobians(toy, state; output_path=path)
        @test comparison.passed
        @test comparison.ordinary_jacobian == comparison.optimized_jacobian == [4.0 3.0; 5.0 2.0]
        @test deserialize(path).completed_columns == 2
        @test deserialize(path).relative_errors == [0.0, 0.0]
        @test state == [2.0, 5.0]
        faulty = values -> begin
            result = toy(values)
            if eltype(values) <: ForwardDiff.Dual{OrdinarySourceDirection}
                result[1] = ForwardDiff.Dual{OrdinarySourceDirection}(
                    ForwardDiff.value(result[1]), (ForwardDiff.partials(result[1])[1] +
                    ForwardDiff.partials(values[1])[1],))
            end
            result
        end
        failed = compare_ad_jacobians(faulty, state)
        @test !failed.passed
        @test failed.absolute_errors == [1.0, 0.0]
        @test failed.completed_columns == 2
        zeros_result = compare_ad_jacobians(values -> [values[1]^2, zero(values[1])], state)
        @test zeros_result.passed && zeros_result.relative_errors[2] == 0.0
        nonfinite = compare_ad_jacobians(values -> [NaN*values[1], values[2]], state)
        @test !nonfinite.passed
        fd = check_finite_differences(toy, state, comparison.optimized_jacobian, values -> true;
            steps=[1e-3, 1e-4, 1e-5], output_path=joinpath(directory, "fd.jls"))
        @test fd.passed && fd.completed_directions == 2
        @test deserialize(joinpath(directory, "fd.jls")).passed
        accepted_path = joinpath(directory, "accepted.jls")
        callback = (iteration, state, mismatch) -> save_verification(accepted_path,
            (; iteration, state, mismatch))
        @test_throws ErrorException newton_solve(values -> [values[1]^2-2], [1.0];
            max_iterations=1, tolerance=1e-12, accepted_state_callback=callback,
            jacobian_backend=:forwarddiff, jacobian_chunk_size=1)
        saved = deserialize(accepted_path)
        @test saved.iteration == 1 && saved.state == [1.5] && saved.mismatch == [0.25]
    end
end

@testset "Full AD matrices through short nonlinear matching" begin
    cloud_m = [1, 2]
    components = [fourier_component(q, [l], cloud_m) for (l, q) in
        [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]]
    grid = angular_grid(8, 9)
    polar_values = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar_values, grid)
    frequencies = ComplexF64[0.405+1e-7im, 0.415+2e-7im]
    amplitudes = [1e-3, 2e-4]
    indices = two_phase_coefficient_indices(components)
    problem = TwoPhaseProblem(; spin=0.7, alpha=0.42, components, cloud_m,
        couplings, polar_values, grid, phase1_nodes=grid.azimuth_nodes,
        phase2_nodes=grid.azimuth_nodes, potential=:cosine, coefficient_scale=2e-4,
        free_outer_coefficient_indices=indices.free_outer_coefficient_indices,
        initial_bound=BitVector([real(sum(component.q .* frequencies)) < 0.42 for component in components]),
        integration=(; inner=4.8, match=5.0, outer=5.2, step=0.1, horizon_order=3, outer_order=7))
    profile = TwoPhaseProfile(frequencies,
        ComplexF64[0.02+0.01im, 0.01-0.003im, 0.002im, -0.001+0.002im],
        ComplexF64[amplitudes..., 1e-5im, -2e-5])
    state = pack_two_phase_profile(profile, problem)
    residual = values -> two_phase_real_matching_residual(values, amplitudes, problem)
    comparison = compare_ad_jacobians(residual, state)
    @test comparison.passed
    @test maximum(comparison.relative_errors) < 1e-10
    @test size(comparison.ordinary_jacobian) == (16, 16)
end
