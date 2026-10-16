using LinearAlgebra: cond, norm

"""
    two_cloud_balance_rates(amplitudes, charge_derivatives, losses; max_condition=1e10)

TODO: describe E/J row ordering, logarithmic column scaling, and signed losses.
"""
function two_cloud_balance_rates(amplitudes, charge_derivatives, losses;
        max_condition=1e10)
    length(amplitudes) == length(losses) == 2 && size(charge_derivatives) == (2, 2) ||
        throw(DimensionMismatch("require a two-cloud E/J balance system"))
    all(value -> isfinite(value) && value > 0, amplitudes) ||
        throw(ArgumentError("require positive finite amplitudes"))
    all(isfinite, charge_derivatives) && all(isfinite, losses) ||
        throw(ArgumentError("require finite derivatives and losses"))
    isfinite(max_condition) && max_condition >= 1 ||
        throw(ArgumentError("require a finite condition limit of at least one"))

    # Columns now multiply d(log A)/dt. Normalize rows to remove the different
    # magnitudes of E and J before judging whether the two equations are independent.
    log_response = charge_derivatives .* transpose(amplitudes)
    row_scales = vec(maximum(abs.(log_response); dims=2))
    all(row_scales .> 0) || error("a charge has no amplitude response")
    scaled_response = log_response ./ row_scales
    scaled_condition = cond(scaled_response)
    isfinite(scaled_condition) && scaled_condition <= max_condition ||
        error("ill-conditioned E/J balance: scaled condition=$scaled_condition")
    log_rates = -(scaled_response \ (losses ./ row_scales))
    amplitude_rates = amplitudes .* log_rates
    all(isfinite, amplitude_rates) && all(isfinite, log_rates) || error("nonfinite amplitude rates")
    balance_residual = charge_derivatives * amplitude_rates + losses
    return (; amplitude_rates, log_rates, scaled_condition,
        raw_condition=cond(charge_derivatives), balance_residual)
end

"""
    two_phase_evolution_sample(problem, amplitudes, state_guess; kwargs...)

TODO: describe profile correction, response evaluation, and total boundary losses.
"""
function two_phase_evolution_sample(problem, amplitudes, state_guess;
        tolerance=1e-9, max_iterations=12, max_condition=1e10)
    solution = correct_two_phase_profile(problem, amplitudes, state_guess;
        tolerance, max_iterations)
    response = two_phase_profile_response(problem, solution; tolerance)
    observables = response.observables
    charges = [observables.energy, observables.angular_momentum]
    losses = [observables.horizon_energy_flux + observables.outer_energy_flux,
        observables.horizon_angular_momentum_flux + observables.outer_angular_momentum_flux]
    balance = two_cloud_balance_rates(amplitudes, response.observable_derivatives,
        losses; max_condition)
    return (; solution, response, balance, charges, losses)
end

"""
    evolve_two_phase(evaluate_sample, initial_sample; final_time, time_step, kwargs...)

TODO: describe midpoint evolution in log amplitudes, profile prediction, and budgets.
"""
function evolve_two_phase(evaluate_sample, initial_sample;
        final_time, time_step, max_log_step=0.05, callback=nothing)
    isfinite(final_time) && final_time >= 0 ||
        throw(ArgumentError("require finite nonnegative final time"))
    isfinite(time_step) && time_step > 0 ||
        throw(ArgumentError("require positive finite time step"))
    isfinite(max_log_step) && max_log_step > 0 ||
        throw(ArgumentError("require positive finite maximum log-amplitude step"))
    all(value -> isfinite(value) && value > 0, initial_sample.solution.amplitudes) ||
        throw(ArgumentError("require positive finite initial amplitudes"))
    step_count = iszero(final_time) ? 0 : ceil(Int, final_time / time_step)
    times = step_count == 0 ? [0.0] : collect(range(0.0, final_time; length=step_count+1))
    samples = [initial_sample]
    integrated_losses = [zeros(2)]
    balance_errors = [zeros(2)]
    !isnothing(callback) && callback(times[1], initial_sample, integrated_losses[1], balance_errors[1])

    for index in 1:step_count
        current = samples[end]
        amplitudes = current.solution.amplitudes
        log_amplitudes = log.(amplitudes)
        spacing = times[index+1] - times[index]
        maximum(abs.(spacing .* current.balance.log_rates)) <= max_log_step ||
            error("time step is too large for the current amplitude rates")

        # Predict the profile using its implicit response, then independently
        # correct it at the midpoint before evaluating the midpoint rate.
        midpoint_amplitudes = exp.(log_amplitudes + spacing .* current.balance.log_rates / 2)
        midpoint_guess = current.solution.state + current.response.state_derivative *
            (midpoint_amplitudes - amplitudes)
        midpoint = evaluate_sample(midpoint_amplitudes, midpoint_guess)
        maximum(abs.(spacing .* midpoint.balance.log_rates)) <= max_log_step ||
            error("time step is too large for the midpoint amplitude rates")
        next_amplitudes = exp.(log_amplitudes + spacing .* midpoint.balance.log_rates)
        next_guess = midpoint.solution.state + midpoint.response.state_derivative *
            (next_amplitudes - midpoint_amplitudes)
        next_sample = evaluate_sample(next_amplitudes, next_guess)

        # Advance the independent boundary-loss integrals with the same
        # midpoint rule. Their comparison with E(t), J(t) tests global balance.
        cumulative = integrated_losses[end] + spacing .* midpoint.losses
        balance_error = next_sample.charges - initial_sample.charges + cumulative
        push!(samples, next_sample)
        push!(integrated_losses, cumulative)
        push!(balance_errors, balance_error)
        !isnothing(callback) && callback(times[index+1], next_sample, cumulative, balance_error)
    end
    return (; times, samples, integrated_losses, balance_errors)
end
