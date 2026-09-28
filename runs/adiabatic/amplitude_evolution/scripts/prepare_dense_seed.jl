# julia --project=. runs/adiabatic/amplitude_evolution/scripts/prepare_dense_seed.jl
# Build a denser cosine seed with ΔA0 <= min(0.25, 0.05 A0) for branch_energy_derivatives.jl.
# Keeps all 178 saved amplitudes and inserts evenly spaced points in wide intervals.
using Serialization
include(joinpath(@__DIR__, "../../../../src/adiabatic/load.jl"))

const campaign = normpath(joinpath(@__DIR__, ".."))
const source_path = joinpath(campaign, "seeds/cosine_branch_20260924.jls")
const dense_seed_path = joinpath(campaign, "seeds/cosine_dense_20260925.jls")
const maximum_amplitude_step = 0.25
const maximum_fractional_step = 0.05

function prepare_dense_seed()
    isfile(dense_seed_path) && error("dense seed already exists; reuse it or choose a new path")
    source = deserialize(source_path)
    source_amplitudes = [entry.amplitude for entry in source.continuation]

    # Corrected endpoint states and state tangents (job 26337413), plus the
    # energy and flux slopes at the same points (job 26400508).
    parts = [deserialize(joinpath(campaign, "results/cosine_derivatives_part$(index).jls"))
        for index in 0:3]
    entries = Dict(entry.branch_index => entry for part in parts for entry in part.solutions)
    endpoint_rows = Dict(row.branch_index => row for row in
        deserialize(joinpath(campaign, "results/cosine_flux_derivatives.jls")).rows)
    for table in (entries, endpoint_rows)
        [table[index].amplitude for index in eachindex(source_amplitudes)] == source_amplitudes ||
            error("derivative tables do not cover the source amplitudes")
    end

    # Cubic Hermite value and slope from endpoint values and slopes.
    hermite(amplitude, left, right, value, slope) = interpolate_branch_energy(amplitude,
        (; amplitude=left.amplitude, energy=getfield(left, value), energy_derivative=getfield(left, slope)),
        (; amplitude=right.amplitude, energy=getfield(right, value), energy_derivative=getfield(right, slope)))

    continuation, predictions = NamedTuple[], NamedTuple[]
    source_indices = Int[]   # 0 marks an inserted point

    for left_index in eachindex(source_amplitudes)
        left = entries[left_index]
        push!(continuation, (; amplitude=left.amplitude, solution=left.solution))
        push!(source_indices, left_index)
        left_index == length(source_amplitudes) && break

        right = entries[left_index+1]
        spacing = right.amplitude - left.amplitude
        allowed_step = min(maximum_amplitude_step, maximum_fractional_step * left.amplitude)
        insert_count = ceil(Int, spacing / allowed_step) - 1

        for insert_index in 1:insert_count
            fraction = insert_index / (insert_count + 1)
            amplitude = left.amplitude + fraction * spacing

            # Hermite state guess; the derivative driver corrects it before use.
            solution = interpolate_branch_state(amplitude, left, right)
            push!(continuation, (; amplitude, solution))
            push!(source_indices, 0)

            # Predictions from the 178-point table, compared after the run.
            left_row, right_row = endpoint_rows[left_index], endpoint_rows[left_index+1]
            energy = hermite(amplitude, left_row, right_row, :energy, :energy_derivative)
            horizon_flux = hermite(amplitude, left_row, right_row, :horizon_flux, :horizon_flux_derivative).energy
            outer_flux = hermite(amplitude, left_row, right_row, :outer_flux, :outer_flux_derivative).energy
            push!(predictions, (; dense_index=length(continuation), amplitude, fraction,
                left_amplitude=left.amplitude, right_amplitude=right.amplitude,
                energy=energy.energy, energy_derivative=energy.energy_derivative,
                horizon_flux, outer_flux, total_flux=horizon_flux+outer_flux,
                amplitude_rate=-(horizon_flux+outer_flux)/energy.energy_derivative))
        end
    end

    serialize(dense_seed_path, merge(source, (; continuation, rows=NamedTuple[],
        source_indices, predictions, maximum_amplitude_step, maximum_fractional_step,
        predictor=:hermite_state)))

    # Checkpoint: originals retained exactly once, and the spacing rule holds.
    amplitudes = [entry.amplitude for entry in continuation]
    step_ratio = maximum(diff(amplitudes) ./ min.(maximum_amplitude_step,
        maximum_fractional_step .* amplitudes[1:end-1]))
    println("Wrote ", dense_seed_path)
    println("points=", length(amplitudes), " original=", count(>(0), source_indices),
        " inserted=", length(predictions))
    println("originals unchanged: ", amplitudes[source_indices .> 0] == source_amplitudes)
    println("strictly increasing: ", all(diff(amplitudes) .> 0))
    println("max step / allowed step: ", step_ratio)
end

prepare_dense_seed()
