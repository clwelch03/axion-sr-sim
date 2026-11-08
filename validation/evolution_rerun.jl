# Rerun of the long fixed-background evolution (scripts/evolution_fixed_background.jl: a = 0.9, α = 0.3, cosine,
# default mesh, A₀ = (1e-6, 1e-7), t ≤ 1.5e9) with the current integrator defaults (ln A state, physics-based first
# step, energy balance as before), compared with the stored trajectory results/fixed_background_a0.9_alpha0.3_solution.csv
# (711 rhs evaluations, state A, Tsit5 heuristic first step, reltol 1e-6).
# Run from the package root:  julia --project=. validation/evolution_rerun.jl
using AxionCloud, Printf, DelimitedFiles
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
stored = readdlm(joinpath(@__DIR__, "..", "results", "fixed_background_a0.9_alpha0.3_solution.csv"), ',', Float64; skipstart=1)
elapsed = @elapsed result = evolve_fixed_background(problem, [1e-6, 1e-7], (0.0, 1.5e9); saveat=stored[:, 1])
trajectory = amplitude_trajectory(result)
@printf("rerun: %d rhs evaluations in %.0f s (stored run: 711 in 3689 s)\n", length(result.diagnostics), elapsed)
println("     t        A211 stored     A211 rerun    rel.diff     A322 stored     A322 rerun    rel.diff")
for index in 1:4:size(stored, 1)
    @printf("%10.3e  %.6e  %.6e  %+.1e  %.6e  %.6e  %+.1e\n", stored[index, 1], stored[index, 2], trajectory.amplitudes[index, 1],
            trajectory.amplitudes[index, 1] / stored[index, 2] - 1, stored[index, 3], trajectory.amplitudes[index, 2],
            trajectory.amplitudes[index, 2] / stored[index, 3] - 1)
end
worst = maximum(abs.(trajectory.amplitudes ./ stored[:, 2:3] .- 1); dims=1)
@printf("largest relative difference over all saved times: A211 %.1e, A322 %.1e\n", worst...)
