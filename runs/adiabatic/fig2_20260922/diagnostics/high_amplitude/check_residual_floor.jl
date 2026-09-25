# Diagnose the high-amplitude matching floor without changing the branch.
ENV["AXION_SEED_PATH"] = get(ENV, "AXION_SEED_PATH",
    joinpath(@__DIR__, "../../seeds/seed_before_residual_check.jls"))
include(joinpath(@__DIR__, "../../../../../src/adiabatic/drivers/continue_cosine.jl"))

function packed_residual(state, amplitude, integration_options)
    unpacked = unpack_complex(state)
    count = length(modes)
    return pack_complex(nonlinear_matching_residual(unpacked[1],
        amplitude .* unpacked[2:count+1],
        amplitude .* vcat(one(unpacked[1]), unpacked[count+2:end]),
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        integration_options...) ./ amplitude)
end

function probe_matching(previous, match_position)
    amplitude, solution = previous.amplitude, previous.solution
    integration_options = merge(options, (; match=match_position))
    state = pack_complex(vcat(solution.frequency,
        solution.horizon_amplitudes ./ amplitude,
        solution.outer_amplitudes[2:end] ./ amplitude))
    residual = state -> packed_residual(state, amplitude, integration_options)
    baseline = residual(state)
    direct = pack_complex(nonlinear_matching_residual(solution.frequency,
        solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        integration_options...) ./ amplitude)
    repeat_difference = norm(residual(state) - baseline)
    spacing = nextfloat.(state) .- state
    differences = Matrix{Float64}(undef, length(baseline), length(state))
    jacobian = similar(differences)
    for column in eachindex(state)
        adjacent = copy(state)
        adjacent[column] = nextfloat(adjacent[column])
        differences[:, column] = residual(adjacent) - baseline
        displacement = 1e-8 * max(1, abs(state[column]))
        plus, minus = copy(state), copy(state)
        plus[column] += displacement
        minus[column] -= displacement
        jacobian[:, column] = (residual(plus) - residual(minus)) / (2displacement)
    end
    observed = [norm(differences[:, column]) for column in eachindex(state)]
    predicted = [norm(jacobian[:, column]) * spacing[column] for column in eachindex(state)]
    column_norms = [norm(jacobian[:, column]) for column in eachindex(state)]
    equilibrated = jacobian ./ transpose(column_norms)
    result = (; match_position, state, baseline, direct,
        repeat_difference, roundtrip_difference=norm(baseline-direct),
        observed, predicted, jacobian, condition=cond(jacobian),
        equilibrated_condition=cond(equilibrated))
    println("PROBE match=", match_position, " residual=", norm(baseline),
        " repeat_difference=", repeat_difference,
        " roundtrip_difference=", result.roundtrip_difference,
        " max_ulp_effect=", maximum(observed), " max_predicted_ulp_effect=", maximum(predicted),
        " condition=", result.condition, " equilibrated_condition=", result.equilibrated_condition)
    println("COMPONENTS match=", match_position, " baseline=", baseline)
    println("ULP_EFFECTS match=", match_position, " observed=", observed, " predicted=", predicted)
    flush(stdout)
    return result
end

function solve_trial(continuation, target, match_position, derivative_increment;
                     prediction=predict_nonlinear_mode(continuation, target))
    integration_options = merge(options, (; match=match_position))
    history = Float64[]
    progress = (iteration, residual_norm) -> begin
        push!(history, residual_norm)
        println("CORRECT match=", match_position, " difference_step=", derivative_increment,
            " iteration=", iteration, " residual=", residual_norm)
        flush(stdout)
    end
    started = time()
    solution = try
        solve_nonlinear_mode(prediction.frequency, target,
            spin, alpha, modes, couplings, polar, azimuthal, grid;
            horizon_guess=prediction.horizon_amplitudes,
            outer_tail_guess=prediction.outer_amplitudes[2:end],
            tolerance=1e-9, max_iterations=12, difference_step=derivative_increment,
            progress_callback=progress, integration_options...)
    catch exception
        exception isa ErrorException || rethrow()
        println("FAILED match=", match_position, " difference_step=", derivative_increment,
            " reason=", sprint(showerror, exception))
        flush(stdout)
        return (; target, integration_options, derivative_increment, history,
            elapsed=time()-started, solution=nothing, error=sprint(showerror, exception))
    end
    profiles = nonlinear_matching_residual(solution.frequency,
        solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        save_profile=true, integration_options...)
    energy = sum(integrate_energy_branch(branch, solution.frequency,
        spin, alpha, modes, polar, azimuthal, grid)
        for branch in (profiles.from_horizon, profiles.from_outer))
    total_flux = outward_flux(profiles.from_outer.states[:,1], options.outer, solution.frequency) -
        outward_flux(profiles.from_horizon.states[:,1], options.inner, solution.frequency)
    result = (; target, integration_options, derivative_increment, history,
        elapsed=time()-started, solution, energy, total_flux)
    println("SOLVED match=", match_position, " difference_step=", derivative_increment,
        " iterations=", solution.iterations, " residual=", solution.scaled_residual_norm,
        " energy=", energy, " total_flux=", total_flux, " frequency=", solution.frequency,
        " elapsed=", result.elapsed)
    flush(stdout)
    return result
end

function run_diagnosis()
    previous = last(saved.continuation)
    println("SEED A0=", previous.amplitude, " saved_residual=", previous.solution.scaled_residual_norm)
    probes, trials = NamedTuple[], NamedTuple[]
    result_path = joinpath(@__DIR__, "residual_floor_diagnosis.jls")
    for match_position in (5.0, 10.0, 20.0)
        push!(probes, probe_matching(previous, match_position))
        serialize(result_path, (; seed_path, probes, trials))
    end
    for (match_position, derivative_increment) in ((5.0, 1e-8), (5.0, 1e-9), (10.0, 1e-8), (20.0, 1e-8))
        push!(trials, solve_trial(saved.continuation, 10.0, match_position, derivative_increment))
        serialize(result_path, (; seed_path, probes, trials))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_diagnosis()
end
