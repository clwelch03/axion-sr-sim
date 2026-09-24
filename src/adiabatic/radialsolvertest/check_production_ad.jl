using Test

repository = normpath(joinpath(@__DIR__, "../../.."))
benchmark_directory = joinpath(repository, "runs/adiabatic/radialsolvertest")

# Load the production driver with the same starting checkpoint as the benchmark.
# Including it defines its functions without starting continuation.
ENV["AXION_SEED_PATH"] = joinpath(benchmark_directory, "seed_cosine.jls")
ENV["AXION_JACOBIAN"] = "forwarddiff"

include(joinpath(
    repository, "runs/adiabatic/preliminary_20260922/branch/continue_cosine.jl",
))

reference = deserialize(joinpath(benchmark_directory, "shooting_ad.jls"))
target = reference.amplitude

# Reproduce the benchmark's prediction, then correct with the production solver.
prediction = predict_nonlinear_mode_cubic(
    saved.continuation, target; history_spacing=0.025,
)

function report_progress(iteration, residual_norm)
    println("iteration=$iteration residual=$residual_norm")
    flush(stdout)
end

solution = solve_nonlinear_mode(
    prediction.frequency, target, spin, alpha, modes, couplings,
    polar, azimuthal, grid;
    horizon_guess=prediction.horizon_amplitudes,
    outer_tail_guess=prediction.outer_amplitudes[2:end],
    jacobian_backend=:forwarddiff,
    tolerance=1e-9, max_iterations=12,
    progress_callback=report_progress, options...,
)

# Compare the converged observables and mode amplitudes.
row = diagnostics(target, solution)

@testset "Production AD agrees with benchmark" begin
    @test solution.scaled_residual_norm <= 1e-9

    for (quantity, tolerance) in (
        (:energy, 1e-7), (:total_flux, 1e-5),
        (:omega_real, 1e-9), (:omega_imag, 1e-5),
    )
        @test isapprox(
            getproperty(row, quantity),
            getproperty(reference.diagnostics, quantity);
            rtol=tolerance,
        )
    end

    @test isapprox(
        solution.horizon_amplitudes, reference.solution.horizon_amplitudes;
        rtol=1e-6,
    )
    @test isapprox(
        solution.outer_amplitudes, reference.solution.outer_amplitudes;
        rtol=1e-6,
    )
end

println(row)