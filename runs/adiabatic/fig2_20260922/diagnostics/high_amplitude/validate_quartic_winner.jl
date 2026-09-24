# Validate saved-method agreement and subdivide the winner's final step.
using Test
ENV["AXION_POTENTIAL"] = "quartic"
ENV["AXION_SEED_PATH"] = joinpath(@__DIR__, "quartic_methods/cubic_h0025_1e10.jls")
include(joinpath(@__DIR__, "../../scripts/continue_polynomial.jl"))

function winning_correction(continuation, target)
    prediction = predict_nonlinear_mode_cubic(continuation, target; history_spacing=0.025)
    history = Float64[]
    progress = (iteration, residual_norm) -> begin
        push!(history, residual_norm)
        println("A0=", target, " iteration=", iteration, " residual=", residual_norm)
        flush(stdout)
    end
    started = time()
    solution = solve_nonlinear_mode(prediction.frequency, target,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        horizon_guess=prediction.horizon_amplitudes,
        outer_tail_guess=prediction.outer_amplitudes[2:end],
        potential=:quartic, tolerance=1e-9, max_iterations=12,
        difference_step=1e-10, progress_callback=progress, options...)
    row = diagnostics(target, solution)
    println("ACCEPTED ", row, " elapsed=", time()-started)
    flush(stdout)
    return (; entry=(; amplitude=target, solution), row, history, elapsed=time()-started)
end

function compare_solutions(label, first_entry, first_row, second_entry, second_row)
    @test first_entry.amplitude == second_entry.amplitude
    @test first_entry.solution.scaled_residual_norm <= 1e-9
    @test second_entry.solution.scaled_residual_norm <= 1e-9
    @test isapprox(first_row.energy, second_row.energy; rtol=1e-7)
    @test isapprox(first_row.omega_real, second_row.omega_real; rtol=1e-9)
    @test isapprox(first_row.omega_imag, second_row.omega_imag; rtol=1e-5)
    @test isapprox(first_row.total_flux, second_row.total_flux; rtol=1e-5)
    @test isapprox(first_entry.solution.horizon_amplitudes,
        second_entry.solution.horizon_amplitudes; rtol=1e-6)
    @test isapprox(first_entry.solution.outer_amplitudes,
        second_entry.solution.outer_amplitudes; rtol=1e-6)
    profiles = map((first_entry, second_entry)) do entry
        solution = entry.solution
        nonlinear_matching_residual(solution.frequency,
            solution.horizon_amplitudes, solution.outer_amplitudes,
            spin, alpha, modes, couplings, polar, azimuthal, grid;
            potential=:quartic, save_profile=true, options...)
    end
    for profile in profiles
        @test norm(profile.residual)/first_entry.amplitude <= 1e-9
    end
    profile_differences = Float64[]
    for side in (:from_horizon, :from_outer)
        first_states = getproperty(profiles[1], side).states
        second_states = getproperty(profiles[2], side).states
        difference = norm(first_states-second_states)/norm(first_states)
        push!(profile_differences, difference)
        @test difference < 1e-6
    end
    comparison = (; label, amplitude=first_entry.amplitude,
        energy_relative=abs(first_row.energy-second_row.energy)/abs(first_row.energy),
        frequency_real_relative=abs(first_row.omega_real-second_row.omega_real)/abs(first_row.omega_real),
        frequency_imag_relative=abs(first_row.omega_imag-second_row.omega_imag)/abs(first_row.omega_imag),
        flux_relative=abs(first_row.total_flux-second_row.total_flux)/abs(first_row.total_flux),
        profile_differences)
    println("COMPARISON ", comparison)
    flush(stdout)
    return comparison
end

@testset "Quartic winner: saved solution and subdivided continuation" begin
    production = deserialize(joinpath(@__DIR__, "../../results/quartic_branch.jls"))
    @test saved.options == production.options
    @test saved.modes == production.modes
    @test saved.potential == production.potential == "quartic"
    first_new = findfirst(entry -> entry.amplitude > saved.initial_amplitude, saved.continuation)
    reference_index = findfirst(entry -> entry.amplitude == saved.continuation[first_new].amplitude,
        production.continuation)
    isnothing(reference_index) && error("Shared production point missing")
    shared_comparison = compare_solutions("existing_production_point",
        saved.continuation[first_new], saved.rows[first_new],
        production.continuation[reference_index], production.rows[reference_index])
    target = last(saved.continuation).amplitude
    history = saved.continuation[1:end-1]
    midpoint = winning_correction(history, (last(history).amplitude+target)/2)
    result_path = joinpath(@__DIR__, "quartic_winner_validation.jls")
    serialize(result_path, (; shared_comparison, midpoint))
    subdivided = winning_correction(vcat(history, [midpoint.entry]), target)
    endpoint_comparison = compare_solutions("subdivided_final_step",
        last(saved.continuation), last(saved.rows), subdivided.entry, subdivided.row)
    for result in (midpoint, subdivided)
        @test result.entry.solution.iterations <= 12
    end
    serialize(result_path, (; shared_comparison, midpoint, subdivided, endpoint_comparison,
        reference_entry=production.continuation[reference_index],
        reference_row=production.rows[reference_index]))
end
