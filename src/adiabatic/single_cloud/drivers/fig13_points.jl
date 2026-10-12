# Batch entry point for the Python figure driver. One Julia process handles
# both radial overtones, so compilation is paid only once.
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 5 || error("usage: fig13_points.jl SEEDS.csv OUTPUT.csv SPIN STEP DECAY_LENGTHS")
    write_growth_rate_points(ARGS[1], ARGS[2]; spin=parse(Float64, ARGS[3]),
        step=parse(Float64, ARGS[4]), decay_lengths=parse(Float64, ARGS[5]))
end
