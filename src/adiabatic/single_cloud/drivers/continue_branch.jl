# julia --project=. src/adiabatic/single_cloud/drivers/continue_branch.jl potential seed.jls output_prefix
# Continue a finite-amplitude stationary branch for the cosine, quartic, or sextic potential.
using Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../load.jl"))

# Command line: positional potential/seed/output; solver settings from AXION_* variables.
function main(arguments)
    length(arguments) == 3 ||
        error("require potential seed.jls output_prefix")
    potential, seed_path, output_prefix = arguments
    haskey(ENV, "AXION_TARGET_AMPLITUDE") || error("set AXION_TARGET_AMPLITUDE")
    environment_float(name, default) = parse(Float64, get(ENV, name, string(default)))

    mkpath(dirname(abspath(output_prefix)))
    continue_branch(
        seed_path, output_prefix;
        potential=Symbol(potential),
        target_amplitude=parse(Float64, ENV["AXION_TARGET_AMPLITUDE"]),
        runtime_seconds=environment_float("AXION_RUN_SECONDS", 1000),
        predictor=Symbol(get(ENV, "AXION_PREDICTOR", "secant")),
        history_spacing=environment_float("AXION_HISTORY_SPACING", 0.025),
        jacobian_backend=Symbol(get(ENV, "AXION_JACOBIAN", "forwarddiff")),
        difference_step=environment_float("AXION_DIFFERENCE_STEP", 1e-8),
        newton_max_iterations=parse(Int, get(ENV, "AXION_NEWTON_MAX_ITERATIONS", "12")),
        maximum_amplitude_step=environment_float("AXION_MAX_AMPLITUDE_STEP", 0.25),
        maximum_fractional_step=environment_float("AXION_MAX_FRACTIONAL_STEP", 0.05),
        initial_step=haskey(ENV, "AXION_INITIAL_STEP") ?
            parse(Float64, ENV["AXION_INITIAL_STEP"]) : nothing,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
