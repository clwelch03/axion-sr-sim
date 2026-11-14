# Build the reference states from scratch and write them to benchmarks/states/<name>.bin.
#
#   julia --project=. benchmarks/build_states.jl [--states linear,plateau,...]
#
# The stored files are the fixed inputs of the regression references. Running this script
# overwrites them with states that differ from the stored ones at the level of the Newton tolerance
# and sometimes more (the linear state, a cold solve, changes in the 3rd digit of some derivatives
# if the probe tolerance changes), so the references must be rewritten afterwards.
include(joinpath(@__DIR__, "common.jl"))
using .BenchmarkHarness
using .BenchmarkHarness: build_state, save_state, selected_specs
using Printf

names = length(ARGS) == 2 && ARGS[1] == "--states" ? split(ARGS[2], ',') : String[]
for spec in selected_specs(names)
    println("building state '$(spec.name)': $(spec.description)"); flush(stdout)
    elapsed = @elapsed state = build_state(spec; verbose=true)
    save_state(state)
    @printf("  built in %.0f s, scaled |F| = %.2e, wrote %s\n", elapsed, state.solution.residual_norm, BenchmarkHarness.state_path(spec))
end
