# Short spindown evolution on the contained baseline mesh (a = 0.9, α = 0.3, cosine, equal initial occupations as in
# scripts/evolution_fixed_background.jl), f_a/M_Pl = 1e-2, t ≤ 2e7: the default action balance against the energy
# balance with the action-based metric work. On a contained mesh the two should agree to ≈ 1e-7 (NOTES.md, Day 3,
# "Metric work, settled"); the energy balance without w is printed as well (expected ≈ 1e-4 off).
# Run from the package root:  julia --project=. validation/spindown_action_check.jl   (≈ 30 min; Slurm)
using AxionCloud, Printf

background = KerrBackground(0.9)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
target_actions = [3e-6, 3e-6]
probe = solve_radial(problem, [1e-6, 1e-6])
amplitudes = probe.amplitudes .* sqrt.(target_actions ./ phase_actions(probe))
@printf("seed: A = (%.6e, %.6e); mesh: %d nodes, r* ∈ [%.2f, %.2f]\n", amplitudes..., length(mesh.rstar), mesh.rstar[1],
        mesh.rstar[end])
flush(stdout)

span = (0.0, 2e7)
decay_constant_ratio = 1e-2
runs = [("action", (; balance=:action)),
        ("energy+w", (; balance=:energy, include_metric_work=true)),
        ("energy, no w", (; balance=:energy, include_metric_work=false))]
finals = Dict{String, Vector{Float64}}()
for (label, options) in runs
    t = @elapsed result = evolve_with_spindown(problem, amplitudes, span, decay_constant_ratio; saveat=1e7, options...)
    trajectory = amplitude_trajectory(result)
    @printf("%-13s %3d rhs evals in %6.1fs\n", label, length(result.diagnostics), t)
    for index in eachindex(trajectory.times)
        @printf("  t = %.3e  A = (%.10e, %.10e)  b = %.15f  j = %.15f\n", trajectory.times[index],
                trajectory.amplitudes[index, 1], trajectory.amplitudes[index, 2], trajectory.black_hole_mass[index],
                trajectory.black_hole_angular_momentum[index])
    end
    finals[label] = vcat(trajectory.amplitudes[end, :], trajectory.black_hole_mass[end],
                         trajectory.black_hole_angular_momentum[end])
    flush(stdout)
end
for label in ("energy+w", "energy, no w")
    reference = finals["action"]
    @printf("%-13s vs action: ΔA/A = (%.3e, %.3e)  Δ(1−b)/(1−b) = %.3e  Δj/(a−j) = %.3e\n", label,
            (finals[label][1:2] ./ reference[1:2] .- 1)...,
            (finals[label][3] - reference[3]) / (1 - reference[3]),
            (finals[label][4] - reference[4]) / (0.9 - reference[4]))
end
