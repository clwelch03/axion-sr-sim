@testset "Primary amplitude target with accelerating growth" begin
    # A'=A^2 has a known doubling time. Its increasing rate also exercises
    # shortening a step whose local time estimate overshoots the target.
    function evaluate(amplitudes, guess)
        state = amplitudes.^2
        derivatives = [2amplitudes[1] 2amplitudes[2]; 2amplitudes[1] 4amplitudes[2]]
        losses = -derivatives * amplitudes.^2
        charges = [sum(state), state[1]+2state[2]]
        solution = (; state, amplitudes=copy(amplitudes))
        response = (; state_derivative=Matrix(Diagonal(2 .* amplitudes)))
        balance = two_cloud_balance_rates(amplitudes, derivatives, losses)
        return (; solution, response, balance, charges, losses)
    end
    initial = evaluate([0.7, 0.3], zeros(2))
    target = 1.4
    expected_time = 1/(2initial.solution.amplitudes[1])
    errors = Float64[]
    for spacing in (0.02, 0.01)
        callbacks = Float64[]
        result = evolve_two_phase_to_amplitude(evaluate, initial;
            target_amplitude=target, time_step=spacing, max_time=1.0,
            relative_tolerance=1e-9,
            callback=(time, sample, losses, errors) -> push!(callbacks, time))
        @test result.termination == :target_reached
        @test abs(result.samples[end].solution.amplitudes[1]/target-1) <= 1e-9
        @test all(sample -> sample.solution.amplitudes[1] <= target*(1+1e-9), result.samples)
        @test result.times == callbacks
        @test all(diff(result.times) .> 0)
        @test result.balance_errors[end] ≈ result.samples[end].charges-initial.charges+
            result.integrated_losses[end]
        push!(errors, abs(result.times[end]-expected_time))
    end
    @test 3 < errors[1]/errors[2] < 5
    limited = evolve_two_phase_to_amplitude(evaluate, initial;
        target_amplitude=target, time_step=0.02, max_time=0.1)
    @test limited.termination == :time_limit
    @test limited.times[end] == 0.1
    short = evolve_two_phase_to_amplitude(evaluate, initial;
        target_amplitude=target, time_step=0.02, max_time=1.0, max_steps=1)
    @test short.termination == :step_limit
    stalled = merge(initial, (; balance=merge(initial.balance, (; log_rates=[0.0, 0.0]))))
    @test evolve_two_phase_to_amplitude(evaluate, stalled;
        target_amplitude=target, time_step=0.02, max_time=1.0).termination == :growth_stalled
    @test_throws ArgumentError evolve_two_phase_to_amplitude(evaluate, initial;
        target_amplitude=0.1, time_step=0.02, max_time=1.0)
    # A restart begins at the saved endpoint; offsets in time and losses are
    # combined by the campaign driver. The amplitude target remains absolute.
    resumed = evolve_two_phase_to_amplitude(evaluate, limited.samples[end];
        target_amplitude=target, time_step=0.02, max_time=1.0,
        relative_tolerance=1e-9)
    @test resumed.termination == :target_reached
    @test abs(resumed.samples[end].solution.amplitudes[1]/target-1) <= 1e-9
end
