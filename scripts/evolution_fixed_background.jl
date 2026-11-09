# Run from the package root:  julia --project=. scripts/evolution_fixed_background.jl   (≈ 1 h; submit to Slurm)
# Long fixed-background evolution into the nonlinear regime (a = 0.9, α = 0.3, cosine): the baseline trajectory for
# the linear / Omiya / joint comparisons. Domain from `domain_estimate` (contains the bound cloud, integrated outer
# data for propagating components), action balance (the default since 2026-11-09), equal initial occupations.
# The day-2 run (default r* ∈ [−80, 180] mesh, series outer boundary, energy balance) is kept in
# results/fixed_background_a0.9_alpha0.3_{diagnostics,solution}.csv; this script writes *_contained_* files.
using AxionCloud, Printf, DelimitedFiles

background = KerrBackground(0.9)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
# Seed both modes with the same occupation (number of quanta), not the same outer coefficient: per unit A the |322⟩
# cloud holds ~1e5 times more quanta than |211⟩, so the day-2 seed A = (1e-6, 1e-7) started with N₃₂₂ ≈ 1000 N₂₁₁.
# At these amplitudes the cloud is linear, N_i ∝ A_i², so one probe solve gives the amplitudes exactly.
target_actions = [3e-6, 3e-6]
probe = solve_radial(problem, [1e-6, 1e-6])
amplitudes = probe.amplitudes .* sqrt.(target_actions ./ phase_actions(probe))
seeded = solve_radial(problem, amplitudes; guess=probe)
@printf("seed: A = (%.6e, %.6e), N = (%.6e, %.6e)\n", amplitudes..., phase_actions(seeded)...)
span = (0.0, 1.5e9)
prefix = joinpath(@__DIR__, "..", "results", "fixed_background_a0.9_alpha0.3_contained")
@printf("mesh: %d nodes, r* ∈ [%.2f, %.2f]; balance = :action\n", length(mesh.rstar), mesh.rstar[1], mesh.rstar[end])
flush(stdout)

t = @elapsed result = evolve_fixed_background(problem, amplitudes, span; saveat=2.5e7, verbose=true)
println("done: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s"); flush(stdout)

open(prefix * "_diagnostics.csv", "w") do io
    println(io, "time,A1,A2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,E,J,L_E,L_J,N1,N2,L_N1,L_N2")
    for d in result.diagnostics
        row = (d.time, d.amplitudes..., real.(d.mode_frequencies)..., imag.(d.mode_frequencies)..., d.charges.energy,
               d.charges.angular_momentum, d.losses.energy, d.losses.angular_momentum, d.actions..., d.action_losses...)
        println(io, join((@sprintf("%.10e", value) for value in row), ","))
    end
end
open(prefix * "_solution.csv", "w") do io
    println(io, "time,A1,A2")
    trajectory = amplitude_trajectory(result)
    for (index, time) in enumerate(trajectory.times)
        @printf(io, "%.8e,%.8e,%.8e\n", time, trajectory.amplitudes[index, 1], trajectory.amplitudes[index, 2])
    end
end
println("wrote $(prefix)_{diagnostics,solution}.csv"); flush(stdout)
