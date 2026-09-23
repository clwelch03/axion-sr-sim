# Compare derivative increments and continuation spacing at a fixed saved seed.
using Test
ENV["AXION_SEED_PATH"] = get(ENV, "AXION_SEED_PATH",
    joinpath(@__DIR__, "seed_before_speed_check.jls"))
include(joinpath(@__DIR__, "continue_cosine.jl"))

function measured_correction(label, continuation, target; difference_step=1e-8)
    prediction = predict_nonlinear_mode(continuation, target)
    history = Float64[]
    progress = (iteration, residual_norm) -> begin
        push!(history, residual_norm)
        println(label, " iteration=", iteration, " scaled_residual=", residual_norm)
        flush(stdout)
    end
    started = time()
    solution = solve_nonlinear_mode(prediction.frequency, target,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        horizon_guess=prediction.horizon_amplitudes,
        outer_tail_guess=prediction.outer_amplitudes[2:end],
        tolerance=1e-9, max_iterations=nothing, difference_step,
        progress_callback=progress, options...)
    solve_seconds = time() - started
    row = diagnostics(target, solution)
    result = (; label, amplitude=target, difference_step,
        initial_residual=first(history), iterations=solution.iterations,
        residual=solution.scaled_residual_norm, energy=row.energy,
        omega_real=row.omega_real, omega_imag=row.omega_imag,
        total_flux=row.total_flux, solve_seconds)
    println("RESULT ", result)
    flush(stdout)
    open(joinpath(@__DIR__, "continuation_speed_validation.csv"), "a") do io
        println(io, join(values(result), ','))
    end
    return (; amplitude=target, solution), row
end

function check_same_solution(first_entry, first_row, second_entry, second_row)
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
end

@testset "Saved branch: derivative increment and larger continuation steps" begin
    open(joinpath(@__DIR__, "continuation_speed_validation.csv"), "w") do io
        println(io, "label,amplitude,difference_step,initial_residual,iterations,residual,energy,omega_real,omega_imag,total_flux,solve_seconds")
    end
    previous = last(saved.continuation)
    println("Baseline A0=", previous.amplitude, " iterations=", previous.solution.iterations)
    # Replay an already accepted point with the same predictor and finer differences.
    replay, replay_row = measured_correction("replay", saved.continuation[1:end-1], previous.amplitude)
    check_same_solution(replay, replay_row, previous, last(saved.rows))
    # Compare one 25% jump with two 12.5% jumps from precisely the same seed.
    target = 1.25previous.amplitude
    direct, direct_row = measured_correction("direct", saved.continuation, target)
    midpoint, midpoint_row = measured_correction("midpoint", saved.continuation, 1.125previous.amplitude)
    intermediate = vcat(saved.continuation, [midpoint])
    subdivided, subdivided_row = measured_correction("subdivided", intermediate, target)
    check_same_solution(direct, direct_row, subdivided, subdivided_row)
end
