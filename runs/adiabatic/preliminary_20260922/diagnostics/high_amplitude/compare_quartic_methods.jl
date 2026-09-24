# Independent, equally budgeted continuation experiments from one frozen seed.
ENV["AXION_POTENTIAL"] = "quartic"
ENV["AXION_SEED_PATH"] = joinpath(@__DIR__, "../../seeds/quartic_before_bounded_correction.jls")
include(joinpath(@__DIR__, "../../branch/continue_polynomial.jl"))

const methods = [
    (label="secant_1e8", predictor=:secant, increment=1e-8, spacing=0.05, step_fraction=1.0),
    (label="secant_1e10", predictor=:secant, increment=1e-10, spacing=0.05, step_fraction=1.0),
    (label="cubic_h001_1e8", predictor=:cubic, increment=1e-8, spacing=0.01, step_fraction=1.0),
    (label="cubic_h001_1e10", predictor=:cubic, increment=1e-10, spacing=0.01, step_fraction=1.0),
    (label="cubic_h0025_1e10", predictor=:cubic, increment=1e-10, spacing=0.025, step_fraction=1.0),
    (label="cubic_h005_1e10", predictor=:cubic, increment=1e-10, spacing=0.05, step_fraction=1.0),
    (label="secant_small_step", predictor=:secant, increment=1e-8, spacing=0.05, step_fraction=0.25),
]

function compare_method(case_index)
    method = methods[case_index]
    result_directory = joinpath(@__DIR__, "quartic_methods")
    mkpath(result_directory)
    result_prefix = joinpath(result_directory, method.label)
    continuation, rows = copy(saved.continuation), copy(saved.rows)
    attempts = NamedTuple[]
    initial_amplitude = last(continuation).amplitude
    amplitude = initial_amplitude
    amplitude_step = method.step_fraction * (5.601310126487904 - initial_amplitude)
    target_amplitude = 5.7
    started = time()

    function save_result(status)
        elapsed = time()-started
        serialize(result_prefix * ".jls.tmp", (; method, seed_path, spin, alpha,
            modes, options, n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
            potential="quartic", tolerance=1e-9, max_iterations=12,
            target_amplitude, initial_amplitude, elapsed, status,
            continuation, rows, attempts))
        mv(result_prefix * ".jls.tmp", result_prefix * ".jls"; force=true)
        summary = (; label=method.label, amplitude, target_amplitude,
            amplitude_gain=amplitude-initial_amplitude, elapsed,
            gain_per_minute=60(amplitude-initial_amplitude)/max(elapsed, eps()),
            accepted=length(continuation)-length(saved.continuation),
            rejected=count(attempt -> !attempt.success, attempts),
            total_corrections=sum((attempt.iterations for attempt in attempts); init=0),
            energy=last(rows).energy, omega_real=last(rows).omega_real,
            omega_imag=last(rows).omega_imag, total_flux=last(rows).total_flux,
            residual=last(rows).residual, status)
        open(result_prefix * ".csv.tmp", "w") do io
            println(io, join(string.(keys(summary)), ','))
            println(io, join(values(summary), ','))
        end
        mv(result_prefix * ".csv.tmp", result_prefix * ".csv"; force=true)
    end

    println("METHOD ", method, " initial_amplitude=", initial_amplitude)
    flush(stdout)
    save_result("running")
    while time()-started < 600 && amplitude < target_amplitude
        target = min(target_amplitude, amplitude + amplitude_step)
        prediction = method.predictor === :cubic ?
            predict_nonlinear_mode_cubic(continuation, target; history_spacing=method.spacing) :
            predict_nonlinear_mode(continuation, target)
        history = Float64[]
        progress = (iteration, residual_norm) -> begin
            push!(history, residual_norm)
            println("A0=", target, " iteration=", iteration, " residual=", residual_norm)
            flush(stdout)
        end
        attempt_started = time()
        failure = ""
        solution = try
            solve_nonlinear_mode(prediction.frequency, target,
                spin, alpha, modes, couplings, polar, azimuthal, grid;
                horizon_guess=prediction.horizon_amplitudes,
                outer_tail_guess=prediction.outer_amplitudes[2:end],
                potential=:quartic, tolerance=1e-9, max_iterations=12,
                difference_step=method.increment, progress_callback=progress, options...)
        catch exception
            exception isa ErrorException || rethrow()
            failure = sprint(showerror, exception)
            println("REJECTED A0=", target, " ", failure)
            nothing
        end
        push!(attempts, (; target, success=!isnothing(solution),
            iterations=max(0, length(history)-1), history, failure,
            elapsed=time()-attempt_started))
        if isnothing(solution)
            amplitude_step /= 2
        else
            row = diagnostics(target, solution)
            push!(continuation, (; amplitude=target, solution))
            push!(rows, row)
            accepted_step = target - amplitude
            amplitude = target
            growth = solution.iterations <= 4 ? 1.6 : solution.iterations <= 8 ? 1.3 : 1.0
            amplitude_step = min(growth*accepted_step, 0.25amplitude)
            println("ACCEPTED ", row)
        end
        flush(stdout)
        save_result("running")
        if amplitude_step < 1e-5amplitude
            save_result("step_exhausted")
            return
        end
    end
    save_result(amplitude >= target_amplitude ? "target_reached" : "budget_reached")
end

case_index = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0")) + 1
compare_method(case_index)
