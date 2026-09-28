# julia --project=. runs/adiabatic/amplitude_evolution/scripts/merge_derivatives.jl [dense-only]
using Serialization, SHA

const derivative_campaign = normpath(joinpath(@__DIR__, ".."))

function read_derivative_shards(seed_path, prefix; shard_count=4)
    seed = deserialize(seed_path)
    input_sha256 = bytes2hex(sha256(read(seed_path)))
    parts = [deserialize(prefix * "_part$(index).jls") for index in 0:shard_count-1]
    reference = first(parts).metadata
    for (index, part) in enumerate(parts)
        metadata = part.metadata
        metadata.input_sha256 == input_sha256 || error("seed hash mismatch: $prefix")
        metadata.shard_index == index - 1 && metadata.shard_count == shard_count ||
            error("shard selection mismatch: $prefix")
        metadata == merge(reference, (; shard_index=index-1)) || error("shard settings differ")
        metadata.flux_derivatives || error("missing flux derivatives")
        for key in (:spin, :alpha, :modes, :n_polar, :n_azimuth)
            getproperty(metadata, key) == getproperty(seed, key) || error("seed settings differ: $key")
        end
        metadata.options == merge(seed.options,
            (; potential=Symbol(seed.potential), fundamental_m=seed.modes[1][2])) ||
            error("seed options differ")
        expected = collect(index:shard_count:length(seed.continuation))
        [row.branch_index for row in part.rows] == expected || error("incomplete/duplicate shard rows")
        [entry.branch_index for entry in part.solutions] == expected || error("incomplete shard states")
        for (row, entry) in zip(part.rows, part.solutions)
            row.amplitude == entry.amplitude == seed.continuation[row.branch_index].amplitude ||
                error("amplitude mismatch")
            all(isfinite, values(row)) || error("nonfinite observable")
            0 <= row.residual <= metadata.tolerance || error("unconverged row")
            row.energy_derivative > 0 || error("energy turning point requires separate treatment")
            all(isfinite, entry.state_derivative) || error("nonfinite state tangent")
            all(isfinite, vcat(entry.solution.frequency, entry.solution.horizon_amplitudes,
                entry.solution.outer_amplitudes)) || error("nonfinite state")
        end
    end
    rows = sort(vcat([part.rows for part in parts]...); by=row -> row.branch_index)
    solutions = sort(vcat([part.solutions for part in parts]...); by=entry -> entry.branch_index)
    return (; seed, rows, solutions, metadata=reference, seed_path, input_sha256)
end

function merge_derivatives(; dense_only=false)
    dense = read_derivative_shards(joinpath(derivative_campaign, "seeds/cosine_dense_20260925.jls"),
        joinpath(derivative_campaign, "results/cosine_dense_derivatives"))
    inputs = [dense]
    if !dense_only
        extension = read_derivative_shards(joinpath(derivative_campaign, "seeds/cosine_extension_20260928.jls"),
            joinpath(derivative_campaign, "results/cosine_extension_derivatives"))
        for key in (:spin, :alpha, :modes, :options, :n_polar, :n_azimuth, :tolerance, :flux_derivatives)
            getproperty(dense.metadata, key) == getproperty(extension.metadata, key) ||
                error("dense and extension settings differ: $key")
        end
        extension.seed.provenance.dense_sha256 == dense.input_sha256 || error("extension uses a different dense seed")
        push!(inputs, extension)
    end
    rows, solutions, provenance = NamedTuple[], NamedTuple[], NamedTuple[]
    for (input_index, input) in enumerate(inputs), (row, entry) in zip(input.rows, input.solutions)
        branch_index = length(rows) + 1
        push!(rows, merge(row, (; branch_index)))
        push!(solutions, merge(entry, (; branch_index)))
        push!(provenance, (; branch_index, input_index, source_branch_index=row.branch_index))
    end
    amplitudes = [row.amplitude for row in rows]
    all(diff(amplitudes) .> 0) || error("overlapping or unordered inputs")
    all(diff([row.energy for row in rows]) .> 0) || error("energy is not monotonic")
    step_ratio = maximum(diff(amplitudes) ./ min.(0.25, 0.05 .* amplitudes[1:end-1]))
    step_ratio <= 1 + 1e-12 || error("dense spacing violated")
    metadata = (; source_metadata=[input.metadata for input in inputs],
        seed_paths=[input.seed_path for input in inputs],
        input_sha256=[input.input_sha256 for input in inputs])
    prefix = joinpath(derivative_campaign, "results/" * (dense_only ? "cosine_dense_derivatives" : "cosine_evolution_derivatives"))
    serialize(prefix * ".jls.tmp", (; metadata, rows, solutions, provenance))
    mv(prefix * ".jls.tmp", prefix * ".jls"; force=true)
    open(prefix * ".csv.tmp", "w") do stream
        println(stream, join(string.(keys(first(rows))), ','))
        for row in rows
            println(stream, join(values(row), ','))
        end
    end
    mv(prefix * ".csv.tmp", prefix * ".csv"; force=true)
    println("Merged ", length(rows), " points, A0=", first(amplitudes), " to ", last(amplitudes),
        "; max residual=", maximum(row.residual for row in rows), "; max step/allowed=", step_ratio)
    return (; metadata, rows, solutions, provenance)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    isempty(ARGS) || ARGS == ["dense-only"] || error("optional argument: dense-only")
    merge_derivatives(; dense_only=!isempty(ARGS))
end
