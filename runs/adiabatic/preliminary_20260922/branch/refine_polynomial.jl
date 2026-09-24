# Fill amplitude gaps using bracketed guesses and independent nonlinear solves.
include(joinpath(@__DIR__, "continue_polynomial.jl"))

function refine_branch()
    spacing = parse(Float64, get(ENV, "AXION_SAMPLE_SPACING", "0.1"))
    runtime = parse(Float64, get(ENV, "AXION_RUN_SECONDS", "13800"))
    isfinite(spacing) && spacing > 0 || error("Require positive finite spacing")
    isfinite(runtime) && runtime > 0 || error("Require positive finite runtime")
    refinement_output = joinpath(@__DIR__, "$(potential_kind)_refined")
    continuation = copy(saved.continuation)
    rows = copy(saved.rows)
    if isfile(refinement_output * ".jls")
        previous = deserialize(refinement_output * ".jls")
        for key in (:potential, :options, :modes, :spin, :alpha, :n_polar, :n_azimuth)
            getproperty(previous, key) == getproperty(saved, key) ||
                error("Refinement checkpoint mismatch: $key")
        end
        amplitudes = Set(entry.amplitude for entry in continuation)
        for (entry, row) in zip(previous.continuation, previous.rows)
            entry.amplitude in amplitudes && continue
            push!(continuation, entry)
            push!(rows, row)
        end
    end
    ordering = sortperm([entry.amplitude for entry in continuation])
    continuation, rows = continuation[ordering], rows[ordering]
    failed_targets = Set{Float64}()
    started = time()
    while time() - started < runtime
        gaps = diff([entry.amplitude for entry in continuation])
        candidates = filter(eachindex(gaps)) do index
            target = (continuation[index].amplitude + continuation[index+1].amplitude)/2
            gaps[index] > spacing && !(target in failed_targets)
        end
        isempty(candidates) && break
        index = candidates[argmax(gaps[candidates])]
        left, right = continuation[index:index+1]
        target = (left.amplitude + right.amplitude)/2
        # Interpolate only the initial guess; save the corrected solution.
        frequency = (left.solution.frequency + right.solution.frequency)/2
        horizon = (left.solution.horizon_amplitudes + right.solution.horizon_amplitudes)/2
        outer = (left.solution.outer_amplitudes + right.solution.outer_amplitudes)/2
        outer[1] = target
        println("FILL A0=$target gap=$(gaps[index]) potential=$potential_kind")
        flush(stdout)
        progress = (iteration, residual) -> begin
            println("  iteration=$iteration residual=$residual")
            flush(stdout)
        end
        solution = try
            solve_nonlinear_mode(frequency, target, spin, alpha, modes,
                couplings, polar, azimuthal, grid;
                horizon_guess=horizon, outer_tail_guess=outer[2:end],
                potential=potential_kind, tolerance=1e-9,
                max_iterations=newton_max_iterations, difference_step, jacobian_backend,
                progress_callback=progress, options...)
        catch exception
            exception isa ErrorException || rethrow()
            println("FAILED A0=$target: ", sprint(showerror, exception))
            push!(failed_targets, target)
            continue
        end
        row = diagnostics(target, solution)
        isfinite(row.residual) && row.residual <= 1e-9 || error("Unconverged refinement")
        insert!(continuation, index+1, (; amplitude=target, solution))
        insert!(rows, index+1, row)
        serialize(refinement_output * ".jls.tmp", (; spin, alpha, modes, options,
            potential=String(potential_kind), n_polar=saved.n_polar,
            n_azimuth=saved.n_azimuth, tolerance=1e-9, difference_step,
            jacobian_backend=String(jacobian_backend),
            newton_max_iterations, spacing, continuation, rows))
        mv(refinement_output * ".jls.tmp", refinement_output * ".jls"; force=true)
        open(refinement_output * ".csv.tmp", "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for saved_row in rows
                println(io, join(values(saved_row), ','))
            end
        end
        mv(refinement_output * ".csv.tmp", refinement_output * ".csv"; force=true)
        println("SAVED ", row)
        flush(stdout)
    end
    println("Finished: points=$(length(rows)), failed_targets=$(collect(failed_targets)), ",
        "largest_gap=$(maximum(diff([entry.amplitude for entry in continuation])))")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    refine_branch()
end
