# Fill amplitude gaps using bracketed guesses and independent nonlinear solves.
include(joinpath(@__DIR__, "../../../../src/adiabatic/single_cloud/load.jl"))

function refine_branch(;
    potential_kind=Symbol(get(ENV, "AXION_POTENTIAL", "quartic")),
    seed_path=joinpath(@__DIR__, "../results/$(potential_kind)_branch.jls"),
    refinement_output=joinpath(@__DIR__, "../results/$(potential_kind)_refined"),
)
    potential_kind in (:quartic, :sextic) || error("require quartic or sextic")
    setup = load_branch_setup(seed_path; potential=potential_kind)
    (; saved, spin, alpha, modes, couplings, polar, azimuthal, grid, options) = setup
    Symbol(saved.potential) == potential_kind || error("seed potential differs")
    jacobian_backend = Symbol(get(ENV, "AXION_JACOBIAN", "finite"))
    difference_step = parse(Float64, get(ENV, "AXION_DIFFERENCE_STEP", "1e-8"))
    newton_max_iterations = parse(Int, get(ENV, "AXION_NEWTON_MAX_ITERATIONS", "12"))
    isfinite(difference_step) && difference_step > 0 || error("require positive finite difference step")
    newton_max_iterations > 0 || error("require positive correction limit")
    spacing = parse(Float64, get(ENV, "AXION_SAMPLE_SPACING", "0.1"))
    runtime = parse(Float64, get(ENV, "AXION_RUN_SECONDS", "13800"))
    isfinite(spacing) && spacing > 0 || error("Require positive finite spacing")
    isfinite(runtime) && runtime > 0 || error("Require positive finite runtime")
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
    mkpath(dirname(refinement_output))
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
        row = branch_diagnostics(setup, target, solution)
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
    return (; continuation, rows)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    refine_branch()
end
