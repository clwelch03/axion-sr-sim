using LinearAlgebra: cond, norm

"""
    two_cloud_balance_rates(amplitudes, charge_derivatives, losses; max_condition=1e10)

Calculate the two cloud amplitude rates from balancing energy and angular momentum.

The charges Q = (E, J) depend on both clouds' amplitudes. Writing `D = (dE/dA_1, dE/dA_2; dJ/dA_1, dJ/dA_2)`
and signed boundary losses `L = (L_E, L_J)`, we have `dQ/dt = D dot{A} = -L` by the chain rule.

Solve for logarithmic rates g_i = d(log A_i)/dt, so that D * Diagonal(amplitudes) * g = -L.
Normalize each row and its corresponding boundary loss by that row's largest absolute coefficient.
This makes the condition number diagnostic independent of the units used for E and J.

Return `log_rates = dot{A} ./ A`, `amplitude_rates = dot{A}`, the scaled
condition number of the log_response matrix, and the residual D dot{A} + L.
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

    # Coefficient matrix for solving for g = d(log A_i)/dt
    # log_response = [D_11 A_1, D_12 A_2; D_21 A_1, D_22 A_2]
    log_response = charge_derivatives .* transpose(amplitudes)

    # Scale log_response such that the condition number check is insensitive to units
    row_scales = vec(maximum(abs.(log_response); dims=2))
    all(row_scales .> 0) || error("a charge has no amplitude response")
    scaled_response = log_response ./ row_scales
    scaled_condition = cond(scaled_response)
    isfinite(scaled_condition) && scaled_condition <= max_condition ||
        error("ill-conditioned E/J balance: scaled condition=$scaled_condition")
    
    # Solve for dot{A}_i/A_i, and then for dot{A}_i itself
    log_rates = -(scaled_response \ (losses ./ row_scales))
    amplitude_rates = amplitudes .* log_rates

    all(isfinite, amplitude_rates) && all(isfinite, log_rates) || error("nonfinite amplitude rates")
    balance_residual = charge_derivatives * amplitude_rates + losses
    return (; log_rates, amplitude_rates, scaled_condition, balance_residual)
end

"""
    evolve_two_phase_to_amplitude(evaluate_sample, initial_sample;
        target_amplitude, time_step, max_time, kwargs...)

TODO: describe the primary-amplitude stopping condition, retry rule, and termination codes.
"""
function evolve_two_phase_to_amplitude(evaluate_sample, initial_sample;
        target_amplitude, time_step, max_time, relative_tolerance=1e-7,
        max_steps=10000, max_log_step=0.05, callback=nothing)
    isfinite(target_amplitude) && target_amplitude > 0 ||
        throw(ArgumentError("require a positive finite target amplitude"))
    isfinite(time_step) && time_step > 0 && isfinite(max_time) && max_time >= 0 ||
        throw(ArgumentError("require a positive time step and nonnegative time limit"))
    isfinite(relative_tolerance) && 0 < relative_tolerance < 1 ||
        throw(ArgumentError("require relative tolerance between zero and one"))
    max_steps isa Integer && max_steps > 0 || throw(ArgumentError("require a positive step limit"))
    isfinite(max_log_step) && max_log_step > 0 || throw(ArgumentError("require a positive log-step limit"))
    all(value -> isfinite(value) && value > 0, initial_sample.solution.amplitudes) ||
        throw(ArgumentError("require positive finite initial amplitudes"))
    initial_sample.solution.amplitudes[1] <= target_amplitude * (1+relative_tolerance) ||
        throw(ArgumentError("initial primary amplitude already exceeds the target"))

    times = [0.0]
    samples = [initial_sample]
    integrated_losses = [zeros(2)]
    balance_errors = [zeros(2)]
    !isnothing(callback) && callback(0.0, initial_sample, integrated_losses[1], balance_errors[1])
    termination = :step_limit
    while true
        current = samples[end]
        remaining_growth = log(target_amplitude / current.solution.amplitudes[1])
        if remaining_growth <= log1p(relative_tolerance)
            termination = :target_reached
            break
        elseif times[end] >= max_time
            termination = :time_limit
            break
        elseif length(times)-1 >= max_steps
            break
        elseif current.balance.log_rates[1] <= 0
            termination = :growth_stalled
            break
        end

        # The local rate estimates the remaining time, but does not decide
        # whether the target was reached. Every accepted endpoint is corrected.
        spacing = min(time_step, max_time-times[end],
            remaining_growth/current.balance.log_rates[1],
            0.8max_log_step/maximum(abs.(current.balance.log_rates)))
        accepted = false
        for attempt in 1:12
            isfinite(spacing) && times[end]+spacing > times[end] ||
                error("no finite positive step remains before the target")
            trial = evolve_two_phase(evaluate_sample, current;
                final_time=spacing, time_step=spacing, max_log_step)
            next_sample = trial.samples[end]
            actual_growth = log(next_sample.solution.amplitudes[1] / current.solution.amplitudes[1])
            if next_sample.solution.amplitudes[1] > target_amplitude * (1+relative_tolerance)
                # If growth accelerated, shorten and recompute the entire step.
                # Interpolating the endpoint would leave an unmatched profile.
                spacing *= clamp(remaining_growth/actual_growth, 0.1, 0.95)
                continue
            end
            cumulative = integrated_losses[end] + trial.integrated_losses[end]
            balance_error = next_sample.charges-initial_sample.charges+cumulative
            push!(times, times[end]+spacing)
            push!(samples, next_sample)
            push!(integrated_losses, cumulative)
            push!(balance_errors, balance_error)
            !isnothing(callback) && callback(times[end], next_sample, cumulative, balance_error)
            accepted = true
            break
        end
        accepted || error("could not resolve the target-amplitude crossing")
    end
    return (; times, samples, integrated_losses, balance_errors, termination)
end

"""
    two_phase_evolution_sample(problem, amplitudes, state_guess; kwargs...)

TODO: describe profile correction, response evaluation, and total boundary losses.
"""
function two_phase_evolution_sample(problem::TwoPhaseProblem, amplitudes, state_guess;
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
