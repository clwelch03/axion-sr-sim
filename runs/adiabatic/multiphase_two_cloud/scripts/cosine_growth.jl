# Joint cosine control: identical seeds, support and time interval to the chiral run.
include(joinpath(@__DIR__, "chiral_growth.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 3 || error("usage: cosine_growth.jl SEED_JLS REFERENCE_TOML OUTPUT_PREFIX")
    run_joint_growth(abspath(ARGS[1]), abspath(ARGS[3]);
        potential=:cosine, reference_path=abspath(ARGS[2]))
end
