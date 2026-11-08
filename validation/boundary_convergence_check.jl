# Are the final r*_max = 400 solutions of validation/boundary_endpoint_table.jl converged? The calibration
# logs non-convergence warnings there (|F| stalling at 1e-7 – 1e-4 in intermediate secant steps, the
# normalization floor of NOTES day 2: |211⟩ is ≈ 1e-8 of its peak near r ≈ 450). Re-solve at the final
# amplitudes from a fresh hydrogenic seed and from the calibrated solution, and report |F| and iterations.
# Also times the outer boundary block of the Jacobian (ForwardDiff through the integration).
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=400.0, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step=3.7)
for outer_data in (:integrated, :series)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, outer_data)
    calibrated = amplitudes_for_energy(problem, [650.0, 650.0])
    @printf("%-10s calibrated: converged = %s  |F| = %.2e  iterations = %d  E = %s  ω = %s\n", outer_data, calibrated.converged,
            calibrated.residual_norm, calibrated.iterations, fundamental_energies(calibrated), mode_frequencies(calibrated))
    again = solve_radial(problem, calibrated.amplitudes; guess=calibrated, reuse_factorization=false)
    @printf("%-10s re-solve:   converged = %s  |F| = %.2e  iterations = %d  Δω = %s\n", outer_data, again.converged, again.residual_norm,
            again.iterations, mode_frequencies(again) .- mode_frequencies(calibrated))
    layout, unknowns, scale, amplitudes = again.layout, again.unknowns, again.scale, again.amplitudes
    parameter_block = unknowns[AC.parameter_range(layout)]
    block = unknowns[AC.node_range(layout, layout.n_nodes)]
    n_state = length(block)
    jacobian_block() = AC.ForwardDiff.jacobian(vcat(block, parameter_block)) do inputs
        AC.boundary_residual(problem, layout, :outer, inputs[1:n_state], inputs[n_state+1:end], amplitudes, scale)
    end
    jacobian_block()
    @printf("%-10s outer-boundary Jacobian block: %.3f s\n", outer_data, @elapsed jacobian_block())
    flush(stdout)
end
