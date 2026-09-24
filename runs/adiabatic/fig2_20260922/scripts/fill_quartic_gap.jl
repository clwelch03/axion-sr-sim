# Continue through the missing interval at explicitly requested amplitudes.
include(joinpath(@__DIR__, "solve_quartic_target.jl"))

function fill_quartic_gap()
    gap_output = joinpath(@__DIR__, "../results/quartic_gap")
    seed_file = joinpath(@__DIR__, "../diagnostics/high_amplitude/quartic_target_4p4.jls")
    gap = deserialize(isfile(gap_output * ".jls") ? gap_output * ".jls" : seed_file)
    refined = deserialize(joinpath(@__DIR__, "../results/quartic_refined.jls"))

    for checkpoint_data in (gap, refined)
        for key in (:potential, :options, :modes, :spin, :alpha, :n_polar, :n_azimuth)
            getproperty(checkpoint_data, key) == getproperty(saved, key) ||
                error("Checkpoint mismatch: $key")
        end
    end

    # Preserve the earlier predictor history, then append our new points in order.
    entries = Dict(entry.amplitude => entry for entry in refined.continuation)
    for entry in saved.continuation
        entries[entry.amplitude] = entry
    end
    history = sort([entry for entry in values(entries) if entry.amplitude < 4.4];
        by=entry -> entry.amplitude)
    append!(history, gap.continuation)
    continuation, rows = copy(gap.continuation), copy(gap.rows)

    # Atomic, separate outputs allow the general refinement job to keep running.
    function save_gap()
        serialize(gap_output * ".jls.tmp", (;
            spin, alpha, modes, options, n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
            potential="quartic", jacobian_backend="forwarddiff", tolerance=1e-9,
            continuation, rows,
        ))
        mv(gap_output * ".jls.tmp", gap_output * ".jls"; force=true)
        open(gap_output * ".csv.tmp", "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for row in rows
                println(io, join(values(row), ','))
            end
        end
        mv(gap_output * ".csv.tmp", gap_output * ".csv"; force=true)
    end

    save_gap()
    for tenth in 45:54
        target = tenth / 10
        target <= last(history).amplitude && continue
        prediction = predict_nonlinear_mode_cubic(history, target; history_spacing=0.05)
        println("Starting quartic A0=$target from ", last(history).amplitude)
        flush(stdout)
        progress = (iteration, residual) -> begin
            println("  iteration=$iteration residual=$residual")
            flush(stdout)
        end

        solution = solve_nonlinear_mode(
            prediction.frequency, target, spin, alpha, modes, couplings,
            polar, azimuthal, grid;
            horizon_guess=prediction.horizon_amplitudes,
            outer_tail_guess=prediction.outer_amplitudes[2:end],
            potential=:quartic, jacobian_backend=:forwarddiff,
            tolerance=1e-9, max_iterations=20,
            progress_callback=progress, options...,
        )
        row = diagnostics(target, solution)
        isfinite(row.residual) && row.residual <= 1e-9 || error("Unconverged target")
        entry = (; amplitude=target, solution)
        push!(history, entry)
        push!(continuation, entry)
        push!(rows, row)
        save_gap()
        println("SAVED ", row)
        flush(stdout)
    end
    println("Finished quartic gap through A0=", last(history).amplitude)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    fill_quartic_gap()
end
