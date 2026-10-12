# Standalone comparison at the amplitude rejected by the original driver.
using Test
ENV["AXION_SEED_PATH"] = get(ENV, "AXION_SEED_PATH",
    joinpath(@__DIR__, "../../seeds/seed_before_secant.jls"))
include(joinpath(@__DIR__, "../../../../../runs/adiabatic/fig2_20260922/scripts/legacy_continue_cosine.jl"))

@testset "Saved branch: extrapolation and converged correction" begin
    target = 0.2782139804271685
    previous = last(saved.continuation)
    ratio = target / previous.amplitude
    scaled = (
        frequency=previous.solution.frequency,
        horizon_amplitudes=ratio .* previous.solution.horizon_amplitudes,
        outer_amplitudes=ratio .* previous.solution.outer_amplitudes,
    )
    extrapolated = predict_nonlinear_mode(saved.continuation, target)
    results = NamedTuple[]
    for (label, prediction) in (("scaled", scaled), ("secant", extrapolated))
        started = time()
        history = Float64[]
        progress = (iteration, residual_norm) -> begin
            push!(history, residual_norm)
            println(label, " iteration=", iteration, " scaled_residual=", residual_norm)
            flush(stdout)
        end
        solution = solve_nonlinear_mode(prediction.frequency, target,
            spin, alpha, modes, couplings, polar, azimuthal, grid;
            horizon_guess=prediction.horizon_amplitudes,
            outer_tail_guess=prediction.outer_amplitudes[2:end],
            tolerance=1e-9, max_iterations=nothing, progress_callback=progress,
            options...)
        row = diagnostics(target, solution)
        @test solution.scaled_residual_norm <= 1e-9
        push!(results, (; label, initial_residual=first(history),
            iterations=solution.iterations, residual=solution.scaled_residual_norm,
            energy=row.energy, omega_real=row.omega_real, omega_imag=row.omega_imag,
            elapsed=time()-started))
        println("RESULT ", last(results))
        flush(stdout)
    end
    @test isapprox(results[1].energy, results[2].energy; rtol=1e-6)
    @test isapprox(results[1].omega_real, results[2].omega_real; rtol=1e-9)
    @test isapprox(results[1].omega_imag, results[2].omega_imag; rtol=1e-5)
    open(joinpath(@__DIR__, "extrapolation_validation.csv"), "w") do io
        println(io, join(string.(keys(first(results))), ','))
        for result in results
            println(io, join(values(result), ','))
        end
    end
end
