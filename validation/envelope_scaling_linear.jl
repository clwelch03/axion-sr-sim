# Envelope-scaled Newton solve on large domains, linear check.
# Free single-mode |211⟩ at a = 0.99, α = 0.3, A = 1e-6 (outer-coefficient normalization). The eigenfrequency
# does not depend on A or on r*_max once the cloud is contained, so every row must reproduce the r*_max = 400 value.
# Before the scaling, r*_max = 1000 collapsed to R ≈ 0 with ω₂₁₁ = 0.29487 (validation/logs/domain_outer_debug_1000.log).
# Cold starts use the hydrogenic seed; warm starts extend the previous (smaller-domain) solution along its tail.
# Run from the package root:  julia --project=. validation/envelope_scaling_linear.jl
using AxionCloud, Printf

background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1)])
amplitudes = [1e-6]
mesh_for(rstar_max) = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0,
                                 stretch_factor=1.03, max_step=3.7)
problem_for(rstar_max) = RadialProblem(background, FreePotential(), 0.3, mode_set; mesh=mesh_for(rstar_max))

reference = solve_radial(problem_for(400.0), amplitudes)
reference_frequency = mode_frequencies(reference)[1]
reference_peak = maximum(abs, profiles(reference)[:, 1])
@printf("reference r*_max = 400: ω = %.13f %+.7ei  converged = %s  iterations = %d\n",
        real(reference_frequency), imag(reference_frequency), reference.converged, reference.iterations)

function report(label, solution)
    frequency = mode_frequencies(solution)[1]
    profile = abs.(profiles(solution)[:, 1])
    @printf("%-28s nodes = %4d  converged = %-5s iterations = %2d  ω = %.13f %+.7ei  |Δω|/|ω| = %.1e  peak/reference = %.8f  tail/peak = %.1e\n",
            label, length(solution.problem.mesh.rstar), solution.converged, solution.iterations, real(frequency),
            imag(frequency), abs(frequency - reference_frequency) / abs(reference_frequency),
            maximum(profile) / reference_peak, profile[end] / maximum(profile))
    flush(stdout)
end

for rstar_max in (600.0, 1000.0, 2000.0, 4000.0)
    report(@sprintf("cold   r*_max = %.0f", rstar_max), solve_radial(problem_for(rstar_max), amplitudes))
end
