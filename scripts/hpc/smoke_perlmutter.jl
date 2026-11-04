# Verify package loading, a transferred checkpoint, analytic averages, and a solve.
using Test, Serialization
repository = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(repository, "src/adiabatic/load.jl"))
include(joinpath(repository, "src/adiabatic/multiphase/load.jl"))
seed = isempty(ARGS) ? joinpath(repository,
    "runs/adiabatic/multiphase_two_cloud/work/qcd_joint_probe_20261030_v2/inputs/primary.jls") : abspath(only(ARGS))
isfile(seed) || error("Transfer the saved inputs first: $seed")
saved = deserialize(seed)
println("Loaded transferred checkpoint: ", seed, " (", typeof(saved), ")")
include(joinpath(repository, "tests/adiabatic/multiphase/test_averaged_observables.jl"))
include(joinpath(repository, "tests/adiabatic/multiphase/test_nonlinear_matching.jl"))
println("PERLMUTTER SMOKE PASSED")
