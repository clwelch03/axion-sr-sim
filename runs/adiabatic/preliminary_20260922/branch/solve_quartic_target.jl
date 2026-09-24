# Solve the requested point inside the remaining quartic sampling gap.
ENV["AXION_POTENTIAL"] = "quartic"
ENV["AXION_JACOBIAN"] = "forwarddiff"
ENV["AXION_SEED_PATH"] = joinpath(@__DIR__, "quartic_branch.jls")
include(joinpath(@__DIR__, "continue_polynomial.jl"))

function solve_quartic_target()
    target = 4.4

    # Combine accepted histories, keeping production values at shared amplitudes.
    refined = deserialize(joinpath(@__DIR__, "quartic_refined.jls"))
    for key in (:potential, :options, :modes, :spin, :alpha, :n_polar, :n_azimuth)
        getproperty(refined, key) == getproperty(saved, key) ||
            error("Refinement checkpoint mismatch: $key")
    end
    entries = Dict(entry.amplitude => entry for entry in refined.continuation)
    for entry in saved.continuation
        entries[entry.amplitude] = entry
    end
    history = sort([entry for entry in values(entries) if entry.amplitude < target];
        by=entry -> entry.amplitude)
    prediction = predict_nonlinear_mode_cubic(history, target; history_spacing=0.05)

    println("Quartic AD: starting from A0=", last(history).amplitude, " toward ", target)
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

    # Keep this output separate while the general refinement job is running.
    result_prefix = joinpath(@__DIR__, "../diagnostics/high_amplitude/quartic_target_4p4")
    serialize(result_prefix * ".jls", (;
        spin, alpha, modes, options, n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
        potential="quartic", jacobian_backend="forwarddiff", tolerance=1e-9,
        continuation=[(; amplitude=target, solution)], rows=[row],
    ))
    open(result_prefix * ".csv", "w") do io
        println(io, join(string.(keys(row)), ','))
        println(io, join(values(row), ','))
    end
    println("SAVED ", row)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    solve_quartic_target()
end
