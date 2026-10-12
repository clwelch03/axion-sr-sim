# One writer per grid case. The Python driver checks hashes and supplies src/.
using Serialization
length(ARGS) == 4 || error("usage: stationary.jl SOURCE_ROOT TARGETS.csv OUTPUT_DIR SOURCE_HASH")
source_root, targets_path, directory, source_hash = ARGS
include(joinpath(source_root, "adiabatic/drivers/qcd_stationary_point.jl"))
mkpath(directory)
for line in readlines(targets_path)[2:end]
    label, numbers... = split(line, ',')
    alpha, spin, energy, step, outer = parse.(Float64, numbers)
    path = joinpath(directory, label * ".jls")
    if isfile(path)
        saved = deserialize(path)
        saved.source_hash == source_hash || error("source changed")
        saved.row.alpha == alpha && saved.row.spin == spin &&
            saved.row.step == step && saved.row.outer == outer &&
            abs(saved.row.energy / energy - 1) < 1e-6 || error("target changed: $label")
        result = saved
    else
        println("START $label alpha=$alpha spin=$spin E=$energy step=$step outer=$outer")
        flush(stdout)
        result = qcd_stationary_point(alpha, spin, energy;
            step, outer, difference_step=min(1e-6, 1e-3 * alpha^3),
            secondary_tolerance=1e-5 * max(1.0, (0.15 / alpha)^2),
            field_diagnostics=true,
        )
        serialize(path * ".tmp", (; source_hash, result...))
        mv(path * ".tmp", path)
    end
    csv_path = joinpath(directory, label * ".csv")
    open(csv_path * ".tmp", "w") do io
        println(io, join(keys(result.row), ','))
        println(io, join(values(result.row), ','))
    end
    mv(csv_path * ".tmp", csv_path; force=true)
    println("ACCEPTED $label")
    flush(stdout)
end
