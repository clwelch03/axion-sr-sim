using Serialization, Printf
campaign = normpath(joinpath(@__DIR__, ".."))
seed = deserialize("$campaign/seeds/cosine_dense_20260925.jls")
include("merge_derivatives.jl")
rows = read_derivative_shards("$campaign/seeds/cosine_dense_20260925.jls",
    "$campaign/results/cosine_dense_derivatives").rows
println("rows=", length(rows), " complete=", [r.branch_index for r in rows] == collect(1:length(seed.continuation)))
println("amplitudes match seed: ", [r.amplitude for r in rows] == [e.amplitude for e in seed.continuation])
println("max residual=", maximum(r.residual for r in rows), " max iterations=", maximum(r.iterations for r in rows))
println("all finite: ", all(r -> all(isfinite, (r.energy_derivative, r.total_flux, r.amplitude_rate)), rows))
println("energy_derivative>0 everywhere: ", all(r -> r.energy_derivative > 0, rows))
println("total CPU s=", round(sum(r.elapsed_seconds for r in rows)))
# originals: recomputed vs 178-point flux table
old = Dict(r.branch_index => r for r in deserialize("$campaign/results/cosine_flux_derivatives.jls").rows)
orig = [(rows[i], old[s]) for (i, s) in enumerate(seed.source_indices) if s > 0]
rel(a, b) = abs(a - b) / abs(b)
for key in (:energy, :energy_derivative, :total_flux, :horizon_flux_derivative, :outer_flux_derivative)
    @printf("original points, max rel diff in %-24s %.2e\n", key, maximum(rel(getfield(n, key), getfield(o, key)) for (n, o) in orig))
end
# inserted: Hermite predictions from the 178-point table vs fresh solves
errors = Dict(k => Float64[] for k in (:energy_derivative, :total_flux, :amplitude_rate))
worst = Dict{Symbol,Any}()
for p in seed.predictions
    r = rows[p.dense_index]; r.amplitude == p.amplitude || error("mismatch")
    for k in keys(errors)
        e = rel(getfield(p, k), getfield(r, k)); push!(errors[k], e)
        (e >= get(worst, k, (0.0,))[1]) && (worst[k] = (e, p.amplitude, p.left_amplitude, p.right_amplitude))
    end
end
n = length(seed.predictions)
for k in (:energy_derivative, :total_flux, :amplitude_rate)
    @printf("inserted n=%d  %-18s within 0.1%%: %3d  max=%.3e%% at A0=%.5g in [%.5g, %.5g]\n", n, k,
        count(<=(1e-3), errors[k]), 100worst[k][1], worst[k][2], worst[k][3], worst[k][4])
end
# where do rate failures sit?
bad = [p.amplitude for (p, e) in zip(seed.predictions, errors[:amplitude_rate]) if e > 1e-3]
println("rate >0.1% at A0 in: ", isempty(bad) ? "none" : (minimum(bad), maximum(bad), length(bad)))
# dense table self-consistency: predict every odd point from its even neighbours (spacing doubled)
