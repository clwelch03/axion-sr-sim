# Envelope-scaled Newton solve on large domains, nonlinear two-mode check.
# Same setup as validation/boundary_endpoint_table.jl (a = 0.99, α = 0.3, cosine, |211⟩ + |322⟩, E₁ = E₂ = 650,
# amplitudes calibrated per domain, integrated outer data, stretched meshes). Rows r*_max = 300 and 400 reproduce
# that table (logs/boundary_endpoint_table.log); r*_max = 600, 800, 1000 reach past the |322⟩ containment radius
# (≈ 500), which the outer-coefficient normalization could not do before the scaling.
# Each row also counts the solve_radial warnings issued during its calibration.
# Run from the package root:  julia --project=. validation/envelope_scaling_joint.jl
using AxionCloud, Printf, Logging

background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh_for(rstar_max) = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0,
                                 stretch_factor=1.03, max_step=3.7)

"""Logger that counts warnings and forwards everything to the console logger."""
mutable struct CountingLogger <: AbstractLogger
    inner::AbstractLogger
    warnings::Int
end
Logging.min_enabled_level(logger::CountingLogger) = Logging.min_enabled_level(logger.inner)
Logging.shouldlog(logger::CountingLogger, args...) = Logging.shouldlog(logger.inner, args...)
function Logging.handle_message(logger::CountingLogger, level, message, args...; kwargs...)
    level == Logging.Warn && (logger.warnings += 1)
    Logging.handle_message(logger.inner, level, message, args...; kwargs...)
end

for rstar_max in (300.0, 400.0, 600.0, 800.0, 1000.0)
    mesh = mesh_for(rstar_max)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    logger = CountingLogger(current_logger(), 0)
    elapsed = @elapsed solution = with_logger(() -> amplitudes_for_energy(problem, [650.0, 650.0]), logger)
    frequencies, losses = mode_frequencies(solution), boundary_losses(solution)
    fundamental_profiles = abs.(profiles(solution)[:, fundamental_channels(mode_set)])
    tail = fundamental_profiles[end, :] ./ vec(maximum(fundamental_profiles; dims=1))
    @printf("r*_max = %4.0f (%4d nodes, %4.0f s, %d warnings, converged = %s): Im ω₂₁₁ = %+.6e  Im ω₃₂₂ = %+.6e  L_E = %.6e  L_J = %.6e  A = (%.6e, %.6e)  Re ω = (%.8f, %.8f)  fundamental tail/peak = (%.1e, %.1e)\n",
            rstar_max, length(mesh.rstar), elapsed, logger.warnings, solution.converged, imag(frequencies[1]),
            imag(frequencies[2]), losses.energy, losses.angular_momentum, solution.amplitudes[1], solution.amplitudes[2],
            real(frequencies[1]), real(frequencies[2]), tail[1], tail[2])
    flush(stdout)
end
