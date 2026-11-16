# Wall time and allocation of each piece of one spindown right-hand-side evaluation.
#
#   julia --project=. benchmarks/timing.jl [--states linear,plateau] [--repeats N] [--skip-evolution] [--csv PATH]
#
# Every piece is called once as a warm-up (compilation, cache fill) and then `N`
# times (default 3, or BENCH_REPEATS); the minimum wall time is reported with the
# bytes allocated by that call. Thread counts come from the environment, as
# `julia -t N` or JULIA_NUM_THREADS and OPENBLAS_NUM_THREADS, and are recorded in the CSV.
include(joinpath(@__DIR__, "common.jl"))
using .BenchmarkHarness
using .BenchmarkHarness: BenchmarkState, load_state, selected_specs, run_metadata, residual_at, jacobian_at, displaced_unknowns,
                         newton_scaling, factorize, evolution_cache, perturbed_solve, same_point_solve, spindown_rhs_call,
                         evolution_segment
using AxionCloud, LinearAlgebra, Printf

struct Measurement
    state::String
    piece::String
    seconds::Float64
    repeats::Int
    allocated_bytes::Int
    note::String
end

"""Minimum wall time of `repeats` calls of `piece()` after a warm-up call, with the bytes allocated by the fastest call."""
function measure(piece::Function, state_name, label; repeats, note=no_note)
    piece()
    best = nothing
    for _ in 1:repeats
        GC.gc()
        timed = @timed piece()
        if best === nothing || timed.time < best.time
            best = timed
        end
    end
    return Measurement(state_name, label, best.time, repeats, best.bytes, note(best.value))
end

no_note(value) = ""
iteration_note(solution) = @sprintf("%d iterations", solution.iterations)

"""All pieces of one state, in the order they appear in a right-hand-side evaluation."""
function measure_state(state::BenchmarkState; repeats)
    name, solution = state.spec.name, state.solution
    displaced = displaced_unknowns(solution)
    matrix = jacobian_at(state, displaced)
    scaling = newton_scaling(state)
    factorization = factorize(matrix, scaling)
    cache = evolution_cache(state)
    directions = amplitude_response(solution)
    layout_size = length(solution.unknowns)
    pieces = [
        ("residual", () -> residual_at(state, displaced), no_note),
        ("jacobian assembly", () -> jacobian_at(state, displaced), no_note),
        ("factorisation", () -> factorize(matrix, scaling), no_note),
        ("back-substitution", () -> factorization \ ones(layout_size), no_note),
        ("solve_radial, perturbed", () -> perturbed_solve(state), iteration_note),
        ("solve_radial, same point", () -> same_point_solve(state), iteration_note),
        ("cloud_charges", () -> cloud_charges(solution), no_note),
        ("phase_actions", () -> phase_actions(solution), no_note),
        ("boundary_losses", () -> boundary_losses(solution), no_note),
        ("action_losses", () -> action_losses(solution), no_note),
        ("amplitude_response", () -> amplitude_response(solution), no_note),
        ("charge_jacobian", () -> charge_jacobian(solution; directions), no_note),
        ("action_jacobian", () -> action_jacobian(solution; directions), no_note),
        ("background_charge_derivatives", () -> background_charge_derivatives(cache, solution, state.background), no_note),
        ("balance_terms (total)", () -> balance_terms(cache, solution, state.background), no_note),
        ("spindown rhs call", () -> spindown_rhs_call(state, cache), no_note),
    ]
    measurements = Measurement[]
    for (label, piece, note) in pieces
        push!(measurements, measure(piece, name, label; repeats, note))
        report(last(measurements))
    end
    return measurements
end

function report(measurement::Measurement)
    @printf("  %-32s %10.4f s %10.1f MB  %s\n", measurement.piece, measurement.seconds, measurement.allocated_bytes / 1e6, measurement.note)
    flush(stdout)
end

"""One run of the evolution segment: wall time, accepted steps and right-hand-side calls."""
function measure_evolution(state::BenchmarkState; repeats)
    best_time, note, bytes = Inf, "", 0
    evolution_segment(state; time_span=(1.5e9, 1.5e9 * (1 + 1e-6)))   # compile with a trivial span
    for _ in 1:repeats
        GC.gc()
        timed = @timed evolution_segment(state)
        result, rhs_calls = timed.value
        if timed.time < best_time
            best_time, bytes = timed.time, timed.bytes
            note = @sprintf("%d accepted steps, %d rhs calls", length(result.ode_solution.t) - 1, rhs_calls)
        end
    end
    return Measurement(state.spec.name, "evolution 1.5e9 -> 1.9e9", best_time, repeats, bytes, note)
end

function write_csv(path, metadata, measurements)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "commit,source_modified,hostname,cpu_model,julia_version,julia_threads,blas_threads,timestamp,",
                "state,piece,seconds_min,repeats,allocated_bytes,note")
        for m in measurements
            fields = (metadata.commit, metadata.source_modified, metadata.hostname, "\"$(metadata.cpu_model)\"", metadata.julia_version,
                      metadata.julia_threads, metadata.blas_threads, metadata.timestamp, m.state, "\"$(m.piece)\"",
                      @sprintf("%.6f", m.seconds), m.repeats, m.allocated_bytes, "\"$(m.note)\"")
            println(io, join(fields, ","))
        end
    end
    println("wrote ", path)
end

function parse_options(arguments)
    options = Dict{Symbol, Any}(:states => String[], :repeats => parse(Int, get(ENV, "BENCH_REPEATS", "3")),
                                :skip_evolution => false, :csv => nothing)
    remaining = collect(arguments)
    while !isempty(remaining)
        argument = popfirst!(remaining)
        if argument == "--states" && !isempty(remaining)
            options[:states] = split(popfirst!(remaining), ',')
        elseif argument == "--repeats" && !isempty(remaining)
            options[:repeats] = parse(Int, popfirst!(remaining))
        elseif argument == "--csv" && !isempty(remaining)
            options[:csv] = popfirst!(remaining)
        elseif argument == "--skip-evolution"
            options[:skip_evolution] = true
        else
            error("unknown argument '$argument'; see the header of benchmarks/timing.jl")
        end
    end
    return options
end

function main(arguments)
    options = parse_options(arguments)
    metadata = run_metadata()
    println("commit $(metadata.commit)$(metadata.source_modified ? " (src/ modified)" : ""), host $(metadata.hostname), ",
            "$(metadata.cpu_model), julia $(metadata.julia_version), julia threads $(metadata.julia_threads), ",
            "BLAS threads $(metadata.blas_threads), repeats $(options[:repeats])")
    measurements = Measurement[]
    for spec in selected_specs(options[:states])
        state = load_state(spec)
        println("state '$(spec.name)': $(length(state.solution.unknowns)) unknowns")
        append!(measurements, measure_state(state; repeats=options[:repeats]))
        if spec.include_evolution && !options[:skip_evolution]
            push!(measurements, measure_evolution(state; repeats=parse(Int, get(ENV, "BENCH_EVOLUTION_REPEATS", "1"))))
            report(last(measurements))
        end
    end
    csv_path = something(options[:csv], joinpath(@__DIR__, "logs",
                         "timing-$(metadata.commit)-$(metadata.hostname)-j$(metadata.julia_threads)b$(metadata.blas_threads).csv"))
    write_csv(csv_path, metadata, measurements)
end

main(ARGS)
