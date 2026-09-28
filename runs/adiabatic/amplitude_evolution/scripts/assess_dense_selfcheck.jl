using Serialization, Printf
include(joinpath(@__DIR__, "../../../../src/adiabatic/load.jl"))
campaign = normpath(joinpath(@__DIR__, ".."))
include("merge_derivatives.jl")
rows = read_derivative_shards("$campaign/seeds/cosine_dense_20260925.jls",
    "$campaign/results/cosine_dense_derivatives").rows
rel(a, b) = abs(a - b) / abs(b)
hermite(A, l, r, v, s) = interpolate_branch_energy(A, (; amplitude=l.amplitude, energy=getfield(l, v), energy_derivative=getfield(l, s)),
                                                      (; amplitude=r.amplitude, energy=getfield(r, v), energy_derivative=getfield(r, s)))
derivative_errors, rate_errors, amplitudes = Float64[], Float64[], Float64[]
for i in 2:length(rows)-1
    l, m, r = rows[i-1], rows[i], rows[i+1]
    e = hermite(m.amplitude, l, r, :energy, :energy_derivative)
    total = hermite(m.amplitude, l, r, :horizon_flux, :horizon_flux_derivative).energy + hermite(m.amplitude, l, r, :outer_flux, :outer_flux_derivative).energy
    push!(derivative_errors, rel(e.energy_derivative, m.energy_derivative)); push!(rate_errors, rel(-total / e.energy_derivative, m.amplitude_rate)); push!(amplitudes, m.amplitude)
end
@printf("dense leave-one-out (doubled spacing), n=%d\n", length(amplitudes))
@printf("  energy_derivative: max %.3e%% at A0=%.4g; within 0.1%%: %d\n", 100maximum(derivative_errors), amplitudes[argmax(derivative_errors)], count(<=(1e-3), derivative_errors))
@printf("  amplitude_rate:    max %.3e%% at A0=%.4g; within 0.1%%: %d\n", 100maximum(rate_errors), amplitudes[argmax(rate_errors)], count(<=(1e-3), rate_errors))
# cosine branch after job 26434700
branch = deserialize(joinpath(campaign, "../fig2_20260922/results/cosine_branch.jls"))
A = [e.amplitude for e in branch.continuation]; F = [r.total_flux for r in branch.rows]
new = A .> 11.9851915227
k = findfirst(i -> F[i] < 0 && F[i+1] >= 0, 1:length(F)-1)
zero_amplitude = A[k] - F[k] * (A[k+1] - A[k]) / (F[k+1] - F[k])
@printf("cosine: points=%d new=%d last A0=%.6g max residual(new)=%.2e\n", length(A), count(new), A[end], maximum(r.residual for r in branch.rows[new]))
@printf("cosine: max step/allowed (new part)=%.3f\n", maximum(diff(A[findfirst(new)-1:end]) ./ min.(0.25, 0.05 .* A[findfirst(new)-1:end-1])))
@printf("cosine: total_flux sign change between A0=%.6f (%.3e) and %.6f (%.3e); linear zero A0=%.5f\n", A[k], F[k], A[k+1], F[k+1], zero_amplitude)
println("cosine: sign changes of total_flux on branch = ", count(i -> sign(F[i]) != sign(F[i+1]), 1:length(F)-1))
println("cosine: energy monotonic in A0 = ", issorted([r.energy for r in branch.rows]))
