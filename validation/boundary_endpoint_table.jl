# Acceptance test of the integrated outer boundary: the joint a = 0.99, α = 0.3 cosine solution at the
# fixed amplitudes A = (4.209e-2, 1.305e-4) (E₁ ≈ E₂ ≈ 650 on the r* ∈ [−100, 200] mesh, NOTES day 2) on
# stretched meshes with r*_max = 200, 300, 400, with the outer data produced by (i) the truncated series
# (the original behaviour), (ii) inward integration for the propagating components, (iii) inward
# integration for every component; then the series-order sensitivity with the new boundary.
# Each row is reached by continuation: free-field solve at A, then cosine solves at A × (1/4, 1/2, 3/4, 1).
# Run from the package root:  julia --project=. validation/boundary_endpoint_table.jl  (≈ 10 min)
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
amplitudes = [4.209e-2, 1.305e-4]
mesh_for(rstar_max) = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step=3.7)

function joint_solution(problem)
    seed = solve_radial(with_potential(problem, FreePotential()), amplitudes)
    solution = seed
    for fraction in (0.25, 0.5, 0.75, 1.0)
        solution = solve_radial(problem, fraction .* amplitudes; guess=solution)
        solution.converged || @warn "not converged at fraction $fraction (|F| = $(solution.residual_norm))"
    end
    return solution
end

function row(label; rstar_max, outer_boundary, horizon_order=3, outer_order=7)
    mesh = mesh_for(rstar_max)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, outer_boundary, horizon_order, outer_order)
    t = @elapsed solution = joint_solution(problem)
    frequencies = mode_frequencies(solution)
    losses = boundary_losses(solution)
    energies = fundamental_energies(solution)
    residual_time = @elapsed AC.residual(problem, solution.layout, solution.unknowns, amplitudes, solution.scale)
    jacobian_time = @elapsed AC.jacobian(problem, solution.layout, solution.unknowns, amplitudes, solution.scale)
    @printf("%-32s r*_max = %3.0f (%4d nodes, %3.0f s): Im ω₂₁₁ = %+.4e  Im ω₃₂₂ = %+.4e  L_E = %.5e  L_J = %.5e  E = (%.2f, %.2f)  Re ω = (%.6f, %.6f)  [residual %.2f s, Jacobian %.2f s]\n",
            label, rstar_max, length(mesh.rstar), t, imag(frequencies[1]), imag(frequencies[2]), losses.energy, losses.angular_momentum,
            energies[1], energies[2], real(frequencies[1]), real(frequencies[2]), residual_time, jacobian_time)
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
row("integrated, all, r*_max = 600"; rstar_max=600.0, outer_boundary=integrated_all)
row("integrated (J_H, J_B) = (5, 9)"; rstar_max=400.0, outer_boundary=integrated, horizon_order=5, outer_order=9)
row("integrated, all, (5, 9)"; rstar_max=400.0, outer_boundary=integrated_all, horizon_order=5, outer_order=9)
row("integrated, all, margin 100, 96/λ"; rstar_max=400.0, outer_boundary=OuterBoundary(include_bound=true, convergence_margin=100.0, steps_per_wavelength=96))
