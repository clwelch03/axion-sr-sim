using Serialization
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

# The batch script supplies the output campaign directory and frozen seed file.
# An optional prefix lets a restart choose fresh files without editing the driver.
function main(arguments)
    length(arguments) in (2, 3) ||
        error("usage: continue_amplitude.jl CAMPAIGN_DIRECTORY SEED_FILE [OUTPUT_PREFIX]")
    campaign = abspath(arguments[1])
    seed = deserialize(abspath(arguments[2]))
    output_prefix = length(arguments) == 3 ? abspath(arguments[3]) :
        joinpath(campaign, "results/amplitude_step_20261028_v1")
    mkpath(dirname(output_prefix))
    return continue_two_phase_amplitude(seed, output_prefix;
        target_amplitude=1.0896058341192524e-4, fractional_step=0.015625,
        tolerance=1e-5, max_iterations=20)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
