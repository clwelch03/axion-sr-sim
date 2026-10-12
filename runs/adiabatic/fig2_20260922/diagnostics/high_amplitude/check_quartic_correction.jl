# Compare direct and subdivided cubic continuation at the formerly slow target.
using Test
ENV["AXION_POTENTIAL"] = "quartic"
ENV["AXION_SEED_PATH"] = joinpath(@__DIR__, "../../seeds/quartic_before_bounded_correction.jls")
include(joinpath(@__DIR__, "../../../../../runs/adiabatic/fig2_20260922/scripts/legacy_continue_polynomial.jl"))

function corrected_point(continuation, target, label)
    prediction = predict_nonlinear_mode_cubic(continuation, target; history_spacing=0.05)
    history = Float64[]
    progress = (iteration, residual_norm) -> begin
        push!(history, residual_norm)
        println(label, " A0=", target, " iteration=", iteration, " residual=", residual_norm)
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
    result = (; label, entry=(; amplitude=target, solution), row, history, elapsed=time()-started)
    println("RESULT ", label, " ", row, " elapsed=", result.elapsed)
    flush(stdout)
    return result
end

@testset "Quartic: bounded correction and continuation spacing" begin
    target = 5.601310126487904
    results = NamedTuple[]
    result_path = joinpath(@__DIR__, "quartic_correction_validation.jls")
    direct = corrected_point(saved.continuation, target, "direct")
    push!(results, direct)
    serialize(result_path, (; seed_path, results))
    midpoint = corrected_point(saved.continuation,
        (last(saved.continuation).amplitude + target)/2, "midpoint")
    push!(results, midpoint)
    serialize(result_path, (; seed_path, results))
    subdivided = corrected_point(vcat(saved.continuation, [midpoint.entry]), target, "subdivided")
    push!(results, subdivided)
    serialize(result_path, (; seed_path, results))
    for result in results
        @test result.entry.solution.scaled_residual_norm <= 1e-9
        @test result.entry.solution.iterations <= 12
    end
    @test isapprox(direct.row.energy, subdivided.row.energy; rtol=1e-7)
    @test isapprox(direct.row.omega_real, subdivided.row.omega_real; rtol=1e-9)
    @test isapprox(direct.row.omega_imag, subdivided.row.omega_imag; rtol=1e-5)
    @test isapprox(direct.row.total_flux, subdivided.row.total_flux; rtol=1e-5)
    @test isapprox(direct.entry.solution.horizon_amplitudes,
        subdivided.entry.solution.horizon_amplitudes; rtol=1e-6)
    @test isapprox(direct.entry.solution.outer_amplitudes,
        subdivided.entry.solution.outer_amplitudes; rtol=1e-6)
end
