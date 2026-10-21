using SHA
include(joinpath(@__DIR__, "response_and_evolution.jl"))

"""
    run_doubling_trajectory(problem, initial, prefix, target, doubling_time, steps)

TODO: describe the absolute amplitude target, checkpoint offsets, and stop reasons.
"""
function run_doubling_trajectory(problem, initial, prefix, target, doubling_time, steps)
    path = prefix * "-steps$(steps)"
    if isfile(path * ".jls")
        saved = deserialize(path * ".jls")
        saved.termination != :running && return saved
        history = (; saved.times, saved.samples, saved.integrated_losses, saved.balance_errors)
    else
        history = (; times=Float64[], samples=typeof(initial)[],
            integrated_losses=Vector{Float64}[], balance_errors=Vector{Float64}[])
    end
    offset_time = isempty(history.times) ? 0.0 : history.times[end]
    offset_losses = isempty(history.times) ? zeros(2) : history.integrated_losses[end]
    starting_sample = isempty(history.times) ? initial : history.samples[end]
    callback = (time, sample, losses, errors) -> begin
        absolute_time = offset_time + time
        !isempty(history.times) && absolute_time == history.times[end] && return
        cumulative = offset_losses + losses
        push!(history.times, absolute_time)
        push!(history.samples, sample)
        push!(history.integrated_losses, cumulative)
        push!(history.balance_errors, sample.charges-initial.charges+cumulative)
        save_checkpoint(path * ".jls", (; history..., termination=:running))
        write_trajectory(path * ".csv", history)
        print_sample("doubling, nominal steps=$steps, t=$absolute_time", sample)
    end
    evaluate = (amplitudes, guess) -> two_phase_evolution_sample(problem,
        amplitudes, guess; tolerance=1e-9)
    result = evolve_two_phase_to_amplitude(evaluate, starting_sample;
        target_amplitude=target, time_step=doubling_time/steps,
        max_time=max(0.0, 2doubling_time-offset_time), relative_tolerance=1e-9, callback)
    finished = (; history..., termination=result.termination)
    save_checkpoint(path * ".jls", finished)
    println("Nominal steps=$steps: ", result.termination, ", final time=", history.times[end])
    flush(stdout)
    return finished
end

"""
    run_primary_doubling(seed_path, prefix)

TODO: describe the frozen seed, amplitude target, and two exponential references.
"""
function run_primary_doubling(seed_path, prefix)
    steps = parse(Int, get(ENV, "AXION_DOUBLING_STEPS", "64"))
    steps > 0 || error("AXION_DOUBLING_STEPS must be positive")
    mkpath(dirname(prefix))
    mkdir(prefix * ".lock")
    try
        manifest_path = joinpath(@__DIR__, "../../../../provenance.json")
        isfile(manifest_path) || error("run from a frozen source snapshot")
        configuration = (; source_manifest=read(manifest_path, String),
            seed_sha256=bytes2hex(sha256(read(seed_path))), steps,
            target_factor=2.0, target_tolerance=1e-9, time_limit_factor=2.0)
        configuration_path = prefix * "-configuration.jls"
        if isfile(configuration_path)
            deserialize(configuration_path) == configuration ||
                error("changed doubling configuration; use a distinct output prefix")
        else
            isfile(prefix * "-initial.jls") && error("existing initial state lacks configuration")
            save_checkpoint(configuration_path, configuration)
        end

        initial_path = prefix * "-initial.jls"
        if isfile(initial_path)
            problem, initial, linear_frequencies = deserialize(initial_path)
        else
            problem, seed = deserialize(seed_path)
            # Re-correct the saved starting profile with the snapshotted code.
            initial = two_phase_evolution_sample(problem, seed.solution.amplitudes,
                seed.solution.state; tolerance=1e-9)
            frequencies = unpack_complex(initial.solution.state)[1:2]
            # The linear solver fixes the series orders internally and accepts
            # only the radial grid keywords from the joint problem.
            problem.integration.horizon_order == 3 && problem.integration.outer_order == 7 ||
                error("linear comparison requires the default boundary orders 3/7")
            linear_grid = (; problem.integration.step, problem.integration.inner,
                problem.integration.outer, problem.integration.match)
            linear_frequencies = map(1:2) do cloud
                label = cloud == 1 ? (1,0) : (0,1)
                component = only(filter(component -> component.q == label, problem.components))
                mode = solve_linear_mode(frequencies[cloud], problem.spin, problem.alpha,
                    component.l_values, problem.cloud_m[cloud], problem.couplings;
                    linear_grid..., tolerance=1e-9)
                mode.frequency
            end
            save_checkpoint(initial_path, (problem, initial, linear_frequencies))
        end
        print_sample("Doubling initial state", initial)
        println("Linear eigenmode frequencies on the same radial grid: ", linear_frequencies)
        initial.balance.log_rates[1] > 0 || error("initial primary is not growing")
        target = 2initial.solution.amplitudes[1]
        doubling_time = log(2) / initial.balance.log_rates[1]
        println("Target A1=", target, "; estimated doubling time=", doubling_time)
        flush(stdout)

        # For an amplitude, the linear growth exponent is Im(omega), not
        # 2Im(omega). The factor of two applies to quadratic occupations/energies.
        open(prefix * "-linear-reference.csv", "w") do io
            println(io, "cloud,initial_amplitude,linear_frequency_real,linear_growth_rate,initial_log_rate")
            for cloud in 1:2
                println(io, join((cloud, initial.solution.amplitudes[cloud],
                    real(linear_frequencies[cloud]), imag(linear_frequencies[cloud]),
                    initial.balance.log_rates[cloud]), ','))
            end
        end
        trajectory = run_doubling_trajectory(problem, initial, prefix, target, doubling_time, steps)
        trajectory.termination == :target_reached || error("run stopped: $(trajectory.termination)")
        final_time = trajectory.times[end]
        linear_expected = initial.solution.amplitudes .* exp.(imag.(linear_frequencies) .* final_time)
        frozen_expected = initial.solution.amplitudes .* exp.(initial.balance.log_rates .* final_time)
        summary = (; target, nominal_steps=steps, final_time,
            deviation_from_linear=trajectory.samples[end].solution.amplitudes ./ linear_expected .- 1,
            deviation_from_initial_rate=trajectory.samples[end].solution.amplitudes ./ frozen_expected .- 1,
            relative_budget_errors=trajectory.balance_errors[end] ./ initial.charges)
        save_checkpoint(prefix * "-summary.jls", summary)
        println("Doubling summary: ", summary)
        println("COMPLETED: primary amplitude doubled.")
    finally
        rm(prefix * ".lock")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("usage: double_primary.jl SEED_JLS OUTPUT_PREFIX")
    run_primary_doubling(abspath(ARGS[1]), abspath(ARGS[2]))
end
