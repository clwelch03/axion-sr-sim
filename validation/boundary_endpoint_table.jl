# Acceptance test of the integrated outer boundary: the joint a = 0.99, α = 0.3 cosine solution at
# E₁ = E₂ = 650 (amplitudes calibrated per domain with `amplitudes_for_energy`, as in NOTES day 2) on
# stretched meshes with r*_max = 200, 300, 400, with the outer data produced by (i) the truncated series
# (the original behaviour), (ii) inward integration for the propagating components, (iii) inward
# integration for every component; then the series-order sensitivity with the new boundary.
# Run from the package root:  julia --project=. validation/boundary_endpoint_table.jl  (≈ 1 h)
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
mesh_for(rstar_max) = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step=3.7)

function row(label; rstar_max, outer_boundary, horizon_order=3, outer_order=7)
    mesh = mesh_for(rstar_max)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, outer_boundary, horizon_order, outer_order)
    t = @elapsed solution = amplitudes_for_energy(problem, [650.0, 650.0])
    frequencies = mode_frequencies(solution)
    losses = boundary_losses(solution)
    amplitudes = solution.amplitudes
    residual_time = @elapsed AC.residual(problem, solution.layout, solution.unknowns, amplitudes, solution.scale)
    jacobian_time = @elapsed AC.jacobian(problem, solution.layout, solution.unknowns, amplitudes, solution.scale)
    @printf("%-34s r*_max = %3.0f (%4d nodes, %4.0f s): Im ω₂₁₁ = %+.4e  Im ω₃₂₂ = %+.4e  L_E = %.5e  L_J = %.5e  A = (%.4e, %.4e)  Re ω = (%.6f, %.6f)  [residual %.2f s, Jacobian %.2f s]\n",
            label, rstar_max, length(mesh.rstar), t, imag(frequencies[1]), imag(frequencies[2]), losses.energy, losses.angular_momentum,
            amplitudes[1], amplitudes[2], real(frequencies[1]), real(frequencies[2]), residual_time, jacobian_time)
    flush(stdout)
end

series = OuterBoundary(data=:series)
integrated = OuterBoundary()
integrated_all = OuterBoundary(include_bound=true)
for rstar_max in (200.0, 300.0, 400.0)
    row("series (original)"; rstar_max, outer_boundary=series)
    row("integrated, propagating only"; rstar_max, outer_boundary=integrated)
    row("integrated, all components"; rstar_max, outer_boundary=integrated_all)
end
row("integrated (J_H, J_B) = (5, 9)"; rstar_max=400.0, outer_boundary=integrated, horizon_order=5, outer_order=9)
row("integrated, all, (5, 9)"; rstar_max=400.0, outer_boundary=integrated_all, horizon_order=5, outer_order=9)
row("integrated, all, margin 100, 96/λ"; rstar_max=400.0, outer_boundary=OuterBoundary(include_bound=true, convergence_margin=100.0, steps_per_wavelength=96))
