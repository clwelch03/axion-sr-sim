# Compare predictors at fixed Jacobian settings using one frozen cosine seed.
using Test
ENV["AXION_SEED_PATH"] = joinpath(@__DIR__, "../../seeds/cosine_before_predictor_comparison.jls")
ENV["AXION_DIFFERENCE_STEP"] = "1e-8"
include(joinpath(@__DIR__, "../../scripts/continue_cosine.jl"))

function compare_cosine_predictors()
    methods = [(label="secant", spacing=0.0),
        [(label="cubic_$(spacing)", spacing=spacing) for spacing in (0.01, 0.025, 0.05, 0.1)]...]
    predict(history, target, method) = method.spacing == 0 ?
        predict_nonlinear_mode(history, target) :
        predict_nonlinear_mode_cubic(history, target; history_spacing=method.spacing)
    initial_amplitude = last(saved.continuation).amplitude
    result_path = joinpath(@__DIR__, "cosine_predictor_comparison.jls")
    screenings, trials = NamedTuple[], NamedTuple[]
    save() = serialize(result_path, (; seed_path, initial_amplitude, difference_step,
        options, modes, screenings, trials))
    println("Frozen seed A0=", initial_amplitude)
    for increment in (0.05, 0.1), method in methods
        target = initial_amplitude + increment
        prediction = predict(saved.continuation, target, method)
        residual = norm(nonlinear_matching_residual(prediction.frequency,
            prediction.horizon_amplitudes, prediction.outer_amplitudes,
            spin, alpha, modes, couplings, polar, azimuthal, grid; options...))/target
        push!(screenings, (; method, increment, target, residual))
        println("SCREEN ", last(screenings))
        flush(stdout)
        save()
    end

    function correct(history, target, method, label)
        prediction = predict(history, target, method)
        residual_history = Float64[]
        progress = (iteration, residual) -> begin
            push!(residual_history, residual)
            println(label, " iteration=", iteration, " residual=", residual)
            flush(stdout)
        end
        started = time()
        failure = ""
        solution = try
            solve_nonlinear_mode(prediction.frequency, target, spin, alpha, modes,
                couplings, polar, azimuthal, grid;
                horizon_guess=prediction.horizon_amplitudes,
                outer_tail_guess=prediction.outer_amplitudes[2:end],
                tolerance=1e-9, max_iterations=12, difference_step,
                progress_callback=progress, options...)
        catch exception
            exception isa ErrorException || rethrow()
            failure = sprint(showerror, exception)
            nothing
        end
        row = isnothing(solution) ? nothing : diagnostics(target, solution)
        trial = (; label, method, target, solution, row, failure,
            elapsed=time()-started, residual_history)
        push!(trials, trial)
        println("RESULT ", label, " elapsed=", trial.elapsed, " row=", row, " failure=", failure)
        flush(stdout)
        save()
        return trial
    end

    # Correct the secant and two best cubic guesses at the same larger target.
    ranked = sort(filter(row -> row.increment == 0.1 && row.method.spacing > 0 &&
        isfinite(row.residual), screenings); by=row -> row.residual)
    candidates = [first(methods), [row.method for row in first(ranked, min(2, length(ranked)))]...]
    direct_trials = [correct(saved.continuation, initial_amplitude+0.1, method,
        "direct_$(method.label)") for method in candidates]
    successful = filter(trial -> !isnothing(trial.solution), direct_trials)
    if isempty(successful)
        println("No 0.1 jump converged within 12 corrections at difference_step=1e-8.")
        return
    end
    winner = successful[argmin([trial.elapsed for trial in successful])]
    midpoint = correct(saved.continuation, initial_amplitude+0.05, winner.method, "midpoint")
    isnothing(midpoint.solution) && return
    history = vcat(saved.continuation, [(amplitude=midpoint.target, solution=midpoint.solution)])
    subdivided = correct(history, winner.target, winner.method, "subdivided")
    isnothing(subdivided.solution) && return
    @testset "Cosine predictor: direct and subdivided larger step" begin
        for trial in (winner, midpoint, subdivided)
            @test trial.solution.scaled_residual_norm <= 1e-9
        end
        for (key, tolerance) in ((:energy, 1e-7), (:omega_real, 1e-9),
                                 (:omega_imag, 1e-5), (:total_flux, 1e-5))
            @test isapprox(getproperty(winner.row, key), getproperty(subdivided.row, key); rtol=tolerance)
        end
        @test isapprox(winner.solution.horizon_amplitudes,
            subdivided.solution.horizon_amplitudes; rtol=1e-6)
        @test isapprox(winner.solution.outer_amplitudes,
            subdivided.solution.outer_amplitudes; rtol=1e-6)
    end
    println("VALIDATED method=", winner.method, " increment=0.1 difference_step=", difference_step)
end

compare_cosine_predictors()
