# Stationary inputs for the preliminary QCD spin-down campaign. This composes
# the existing Omiya solvers; the primary is nonlinear and the secondary uses
# the phase-averaged primary background. No BH time evolution occurs here.
using Serialization, SHA
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    # Header: label,alpha,spin,energy,step,outer. One writer per output folder.
    length(ARGS) == 2 || error("usage: qcd_stationary_point.jl TARGETS.csv OUTPUT_DIRECTORY")
    targets_path, directory = ARGS
    mkpath(directory)
    # The driver delegates the solve to model libraries. Include their contents
    # in restart provenance so changing the solver invalidates an old result.
    source_root = normpath(joinpath(@__DIR__, "../.."))
    source_files = sort([
        joinpath(directory, filename)
        for (directory, _, filenames) in walkdir(source_root)
        for filename in filenames if endswith(filename, ".jl")
    ])
    source_hash = bytes2hex(sha256(join([
        relpath(path, source_root) * ":" * bytes2hex(sha256(read(path)))
        for path in source_files
    ], "\n")))
    Threads.@threads for line in readlines(targets_path)[2:end]
        label, numbers... = split(line, ',')
        alpha, spin, energy, step, outer = parse.(Float64, numbers)
        path = joinpath(directory, label * ".jls")
        if isfile(path)
            saved = deserialize(path)
            saved.source_hash == source_hash || error("source changed; use a new output folder")
            saved.row.alpha == alpha && saved.row.spin == spin &&
                saved.row.step == step && saved.row.outer == outer &&
                abs(saved.row.energy / energy - 1) < 1e-6 ||
                error("target changed for $label; use a new output folder")
            continue
        end
        println("START $label $numbers")
        flush(stdout)
        result = qcd_stationary_point(alpha, spin, energy; step, outer)

        # Each label has one writer. Publish the binary checkpoint atomically
        # and retain the source hash so incompatible resumes are rejected.
        serialize(path * ".tmp", (; source_hash, result...))
        mv(path * ".tmp", path)
        open(joinpath(directory, label * ".csv"), "w") do io
            println(io, join(keys(result.row), ','))
            println(io, join(values(result.row), ','))
        end
        println("ACCEPTED $label $(result.row)")
        flush(stdout)
    end
end
