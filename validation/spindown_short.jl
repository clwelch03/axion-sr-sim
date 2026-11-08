# Short spindown evolution in the linear regime with f_a/M_Pl = 1e-2: end-to-end check that the hole
# loses (f_a/M_Pl)² × ΔE_cloud and that A₂₁₁ tracks the fixed-background growth.
# Run from the package root:  julia --project=. validation/spindown_short.jl  (≈ 20 min)
using AxionCloud, Printf
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
amplitudes = [1e-6, 1e-7]
t = @elapsed result = evolve_with_spindown(problem, amplitudes, (0.0, 2e7), 1e-2; saveat=1e7, verbose=true)
println("spindown run: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s")
trajectory = amplitude_trajectory(result)
for index in eachindex(trajectory.times)
    @printf("t = %.3e  A = (%.6e, %.6e)  b = %.15f  j = %.15f\n", trajectory.times[index], trajectory.amplitudes[index, 1],
            trajectory.amplitudes[index, 2], trajectory.black_hole_mass[index], trajectory.black_hole_angular_momentum[index])
end
growth = imag(result.diagnostics[1].mode_frequencies[1])
@printf("fixed-background prediction at t = 2e7: A1 = %.6e\n", amplitudes[1] * exp(growth * 2e7))
