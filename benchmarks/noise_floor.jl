# Measure the noise floor of every regression quantity, to set the tolerances in tolerances.toml.
#
#   julia --project=. benchmarks/noise_floor.jl [--states linear,plateau] [--evolution] [--csv PATH]
#
# The baseline is the stored state collected as `regression.jl` does. Each variation changes
# something that should not change the numbers; the relative difference of a quantity from its baseline
# (the measure defined in quantities.jl) is its noise under that variation. Two tiers of quantity exist
# (see README.md): evaluations at the stored unknown vector, and quantities after a re-solve.
#
#   solver_path     the perturbed solve and the evolution use a ten times smaller Newton tolerance (not on the
#                   1355-node mesh, whose residual floor is 5e-10), and the perturbed solve starts from a
#                   solution at A2 x 1.02 instead of from the state. Sets the re-solve tier.
#   roundoff_input  every unknown is multiplied by 1 + 4e-16 x (fixed pattern): how much roundoff-level
#                   differences of the input are amplified. Estimates the effect of reordered floating-point
#                   sums for the evaluation tier.
#   blas_threads    BLAS with 4 threads instead of 1.
#
# A different CPU model cannot be varied in one process: write the reference on one machine and run
# `regression.jl` on another. The class floors printed at the end are the maxima over variations and states.
include(joinpath(@__DIR__, "common.jl"))
using .BenchmarkHarness
using .BenchmarkHarness: collect_quantities, relative_difference, load_state, selected_specs,
                         perturbed_amplitudes, problem_at, run_metadata, solution_from_unknowns, fixed_vector, BenchmarkState
using AxionCloud, LinearAlgebra, Printf

struct Variation
    name::String
    collect::Function   # state -> QuantityTable
end

function variations(evolution::Bool)   # evolution: also run the segment on states that have one
    solver_path = state -> begin
        tighter = state.spec.name == "large_mesh" ? state.spec.run_tolerance : state.spec.run_tolerance / 10
        neighbour = solve_radial(state.problem, state.solution.amplitudes .* [1.0, 1.02]; guess=state.solution, tol=state.spec.run_tolerance)
        collect_quantities(state; newton_tolerance=tighter, evolution, perturbed_guess=neighbour)
    end
    roundoff_input = state -> begin   # unknowns changed by a few ulp, a proxy for reordered floating-point sums
        displaced = state.solution.unknowns .* (1 .+ 4e-16 .* fixed_vector(length(state.solution.unknowns); frequency=0.91))
        solution = solution_from_unknowns(state.problem, state.solution.amplitudes, displaced, state.solution.scale)
        collect_quantities(BenchmarkState(state.spec, state.reference, state.problem, state.background, solution); evolution)
    end
    blas_threads = state -> begin
        previous = BLAS.get_num_threads()
        BLAS.set_num_threads(4)
        try
            collect_quantities(state; evolution)
        finally
            BLAS.set_num_threads(previous)
        end
    end
    return [Variation("solver_path", solver_path), Variation("roundoff_input", roundoff_input), Variation("blas_threads", blas_threads)]
end

function main(arguments)
    state_names, evolution, csv_path = String[], false, nothing
    remaining = collect(arguments)
    while !isempty(remaining)
        argument = popfirst!(remaining)
        if argument == "--states" && !isempty(remaining)
            state_names = split(popfirst!(remaining), ',')
        elseif argument == "--evolution"
            evolution = true
        elseif argument == "--csv" && !isempty(remaining)
            csv_path = popfirst!(remaining)
        else
            error("unknown argument '$argument'; see the header of benchmarks/noise_floor.jl")
        end
    end
    BLAS.set_num_threads(1)
    rows = NamedTuple[]
    for spec in selected_specs(state_names)
        state = load_state(spec)
        with_evolution = evolution && spec.include_evolution
        baseline = collect_quantities(state; evolution=with_evolution)
        for variation in variations(with_evolution)
            variation_table = variation.collect(state)
            for ((name, reference), (_, varied)) in zip(baseline, variation_table)
                difference = relative_difference(varied.values, reference.values, reference.anchor)
                push!(rows, (state=spec.name, variation=variation.name, quantity=name, class=reference.class, difference))
            end
            @printf("state '%s', variation '%s': max difference %.2e\n", spec.name, variation.name,
                    maximum(row.difference for row in rows if row.state == spec.name && row.variation == variation.name))
            flush(stdout)
        end
    end
    println("\nclass                  variation        state         max rel. diff   worst quantity")
    for class in sort(unique(row.class for row in rows)), variation in unique(row.variation for row in rows)
        selected = [row for row in rows if row.class == class && row.variation == variation]
        isempty(selected) && continue
        worst = argmax(row -> isfinite(row.difference) ? row.difference : Inf, selected)
        @printf("%-22s %-16s %-12s %12.3e   %s\n", class, variation, worst.state, worst.difference, worst.quantity)
    end
    println("\nfloor per class (max over variations and states)")
    for class in sort(unique(row.class for row in rows))
        @printf("%-22s %12.3e\n", class, maximum(row.difference for row in rows if row.class == class))
    end
    if csv_path !== nothing
        open(csv_path, "w") do io
            println(io, "state,variation,quantity,class,relative_difference")
            for row in rows
                @printf(io, "%s,%s,%s,%s,%.6e\n", row.state, row.variation, row.quantity, row.class, row.difference)
            end
        end
    end
end

main(ARGS)
