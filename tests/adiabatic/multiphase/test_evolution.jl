using LinearAlgebra: Diagonal

@testset "E/J balance scaling and singularity" begin
    amplitudes = [0.3, 0.02]
    derivatives = [2amplitudes[1] 2amplitudes[2]; 2amplitudes[1] 4amplitudes[2]]
    expected_rates = [0.04, -0.001]
    losses = -derivatives * expected_rates
    result = two_cloud_balance_rates(amplitudes, derivatives, losses)
    @test result.amplitude_rates ≈ expected_rates
    @test result.log_rates ≈ expected_rates ./ amplitudes
    @test norm(result.balance_residual) < 1e-15
    # Changing the units of E or J must not change the solution or its
    # dimensionless condition diagnostic.
    row_units = [1e12, 1e-9]
    scaled = two_cloud_balance_rates(amplitudes, derivatives .* row_units, losses .* row_units)
    @test scaled.amplitude_rates ≈ expected_rates
    @test scaled.scaled_condition ≈ result.scaled_condition
    @test_throws ErrorException two_cloud_balance_rates(amplitudes, [1.0 2; 2 4], losses)
    @test_throws ArgumentError two_cloud_balance_rates([0.0, 1.0], derivatives, losses)
    @test_throws DimensionMismatch two_cloud_balance_rates(amplitudes, zeros(2,3), losses)
end

@testset "Midpoint evolution: analytic amplitudes and independent charge budgets" begin
    # Manufactured dynamics dA_i/dt=A_i^2 has exact solution A_i(0)/(1-t*A_i(0)).
    # Both E=A1^2+A2^2 and J=A1^2+2A2^2 vary, so their flux integrals test
    # a different observable from the amplitudes themselves.
    function evaluate(amplitudes, guess)
        all(isfinite, guess) || error("nonfinite predictor")
        state = amplitudes.^2
        derivatives = [2amplitudes[1] 2amplitudes[2]; 2amplitudes[1] 4amplitudes[2]]
        losses = -derivatives * amplitudes.^2
        charges = [sum(state), state[1]+2state[2]]
        solution = (; state, amplitudes=copy(amplitudes))
        response = (; state_derivative=Matrix(Diagonal(2 .* amplitudes)))
        balance = two_cloud_balance_rates(amplitudes, derivatives, losses)
        return (; solution, response, balance, charges, losses)
    end
    initial_amplitudes = [0.7, 0.3]
    initial = evaluate(initial_amplitudes, initial_amplitudes.^2)
    final_time = 0.2
    exact = initial_amplitudes ./ (1 .- final_time .* initial_amplitudes)
    amplitude_errors, budget_errors = Float64[], Float64[]
    for steps in (4, 8, 16)
        callbacks = Float64[]
        result = evolve_two_phase(evaluate, initial; final_time,
            time_step=final_time/steps, max_log_step=0.1,
            callback=(time, sample, losses, errors) -> push!(callbacks, time))
        push!(amplitude_errors, norm(result.samples[end].solution.amplitudes-exact))
        push!(budget_errors, norm(result.balance_errors[end]))
        @test callbacks == result.times
        @test length(result.samples) == steps+1
        @test all(sample -> all(sample.solution.amplitudes .> 0), result.samples)
        @test result.balance_errors[end] ≈ result.samples[end].charges-initial.charges +
            result.integrated_losses[end]
    end
    @test all(3.5 .< amplitude_errors[1:2] ./ amplitude_errors[2:3] .< 4.5)
    @test all(3.5 .< budget_errors[1:2] ./ budget_errors[2:3] .< 4.5)
    @test initial.solution.amplitudes == initial_amplitudes
    unchanged = evolve_two_phase(evaluate, initial; final_time=0.0, time_step=0.1)
    @test unchanged.samples == [initial]
    @test_throws ErrorException evolve_two_phase(evaluate, initial; final_time=1.0, time_step=1.0)
    @test_throws ArgumentError evolve_two_phase(evaluate, initial; final_time=-1.0, time_step=0.1)
    @test_throws ArgumentError evolve_two_phase(evaluate, initial; final_time=1.0, time_step=0.0)
end
