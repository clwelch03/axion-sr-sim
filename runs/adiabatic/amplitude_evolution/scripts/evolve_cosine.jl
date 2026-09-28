# Run from the repository root with julia --project=.
# Evolve the 372 dense + 54 extension points on their saved spatial grid.
using Serialization, SHA
include(joinpath(@__DIR__, "../../../../src/adiabatic/load.jl"))

function evolve_cosine()
    campaign = normpath(joinpath(@__DIR__, ".."))
    source_path = joinpath(campaign, "results/cosine_evolution_derivatives.jls")
    source_bytes = read(source_path)
    data = deserialize(IOBuffer(source_bytes))
    @assert length(data.rows) == 426
    @assert count(entry -> entry.input_index == 1, data.provenance) == 372
    @assert count(entry -> entry.input_index == 2, data.provenance) == 54

    branch = (;
        amplitudes=[row.amplitude for row in data.rows],
        samples=data.rows,
    )
    initial_amplitude = first(branch.amplitudes)
    final_time = 4e7
    time_step = 250.0
    comparison_time_step = 500.0

    coarse = evolve_amplitude(
        branch, initial_amplitude; final_time, time_step=comparison_time_step)
    fine = evolve_amplitude(branch, initial_amplitude; final_time, time_step)

    # Compare the whole trajectory at shared times, not just the endpoint.
    @assert coarse.times == fine.times[1:2:end]
    relative_difference = maximum(abs.(
        coarse.amplitudes ./ fine.amplitudes[1:2:end] .- 1))
    @assert relative_difference < 1e-5
    @assert all(diff(fine.amplitudes) .>= 0)

    println("Maximum timestep discrepancy: ", 100relative_difference, "%")
    println("Final amplitude: ", last(fine.amplitudes))

    output_prefix = joinpath(campaign, "results/cosine_evolution")
    metadata = (;
        source_path,
        source_sha256=bytes2hex(sha256(source_bytes)),
        initial_amplitude,
        final_time,
        time_step,
        comparison_time_step,
        maximum_relative_difference=relative_difference,
    )
    serialize(output_prefix * ".jls.tmp", (; metadata, trajectory=fine))
    mv(output_prefix * ".jls.tmp", output_prefix * ".jls"; force=true)

    open(output_prefix * ".csv.tmp", "w") do stream
        println(stream, "time_over_M,amplitude,amplitude_rate")
        for (time, amplitude) in zip(fine.times, fine.amplitudes)
            rate = amplitude_rate(amplitude, branch)
            println(stream, time, ',', amplitude, ',', rate)
        end
    end
    mv(output_prefix * ".csv.tmp", output_prefix * ".csv"; force=true)
    println("Saved ", output_prefix, ".{csv,jls}")
    return (; metadata, trajectory=fine)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    evolve_cosine()
end
