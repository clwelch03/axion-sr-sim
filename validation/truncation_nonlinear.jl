# Truncation convergence of the nonlinear rates at a strongly deformed two-mode state.
# a = 0.99, α = 0.3, cosine, |211⟩ + |322⟩ calibrated to E₁ = E₂ = 650 (fundamental energies), on a domain that
# contains the cloud (r*_max = 600, possible since the envelope scaling of the Newton system). Varies the retained
# Fourier order (max_order 3 → 5) and the polar indices per component (n_angular 1 → 2). Production: (3, 1).
# The earlier checks (NOTES day 2: radial step, polar and phase quadrature) did not vary either truncation.
# Run from the package root:  julia --project=. validation/truncation_nonlinear.jl
using AxionCloud, Printf

background = KerrBackground(0.99)
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=600.0, step=0.4, stretch_from=150.0, stretch_factor=1.03,
                  max_step=3.7)
modes = [BoundMode(2, 1, 1), BoundMode(3, 2, 2)]

for (max_order, n_angular) in ((3, 1), (5, 1), (3, 2), (5, 2))
    mode_set = ModeSet(modes; max_order, n_angular)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    elapsed = @elapsed solution = amplitudes_for_energy(problem, [650.0, 650.0])
    frequencies, losses = mode_frequencies(solution), boundary_losses(solution)
    action_rates = action_losses(solution)
    @printf("max_order = %d, n_angular = %d (%2d components, %3d channels, %4.0f s, converged = %s): Im ω₂₁₁ = %+.6e  Im ω₃₂₂ = %+.6e  Re ω = (%.9f, %.9f)  L_E = %.6e  L_J = %.6e  L_N = (%+.6e, %+.6e)  A = (%.6e, %.6e)\n",
            max_order, n_angular, length(mode_set.components), length(mode_set.channel_l), elapsed, solution.converged,
            imag(frequencies[1]), imag(frequencies[2]), real(frequencies[1]), real(frequencies[2]), losses.energy,
            losses.angular_momentum, action_rates[1], action_rates[2], solution.amplitudes[1], solution.amplitudes[2])
    flush(stdout)
end
