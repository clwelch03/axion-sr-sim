# Freeze only the newly continued cosine points for the derivative driver.
using Serialization, SHA

const campaign = normpath(joinpath(@__DIR__, ".."))
const source_path = joinpath(campaign, "../fig2_20260922/results/cosine_branch.jls")
const dense_path = joinpath(campaign, "seeds/cosine_dense_20260925.jls")
const output_path = joinpath(campaign, "seeds/cosine_extension_20260928.jls")

function prepare_extension_seed()
    isfile(output_path) && error("extension seed already exists; reuse it")
    source, dense = deserialize(source_path), deserialize(dense_path)
    for key in (:spin, :alpha, :modes, :options, :n_polar, :n_azimuth)
        getproperty(source, key) == getproperty(dense, key) || error("settings differ: $key")
    end
    Symbol(source.potential) == Symbol(dense.potential) == :cosine || error("require cosine")
    dense_amplitudes = [entry.amplitude for entry in dense.continuation]
    source_amplitudes = [entry.amplitude for entry in source.continuation]
    all(diff(source_amplitudes) .> 0) || error("source amplitudes must increase")
    all(amplitude -> amplitude in dense_amplitudes,
        filter(<=(last(dense_amplitudes)), source_amplitudes)) || error("source prefix differs")
    source_indices = findall(>(last(dense_amplitudes)), source_amplitudes)
    isempty(source_indices) && error("no new points")
    continuation = source.continuation[source_indices]
    amplitudes = vcat(last(dense_amplitudes), [entry.amplitude for entry in continuation])
    step_ratio = maximum(diff(amplitudes) ./ min.(0.25, 0.05 .* amplitudes[1:end-1]))
    step_ratio <= 1 + 1e-12 || error("extension violates dense spacing")
    provenance = (; source_path, source_sha256=bytes2hex(sha256(read(source_path))),
        dense_path, dense_sha256=bytes2hex(sha256(read(dense_path))), source_indices)
    serialize(output_path, merge(source, (; continuation, rows=source.rows[source_indices], provenance)))
    println("Wrote ", output_path, ": ", length(continuation), " points, A0=",
        first(continuation).amplitude, " to ", last(continuation).amplitude,
        "; max step/allowed=", step_ratio)
end

prepare_extension_seed()
