# Reproduce the joint-deformation comparison: |R_{l,q}| of the primary |211⟩ and
# secondary |322⟩ clouds versus r*/M at equal per-cloud energy, solved three ways:
#   Linear  – free field,
#   Single  – each cloud alone with self-interactions (the "Omiya" curves),
#   Joint   – both clouds coupled through the nonlinear source.
# Run from the package root:  julia --project=scripts scripts/profile_comparison.jl
using AxionCloud, CairoMakie, Printf

spin = 0.99
mass_coupling = 0.3
target_energy = 650.0
potential = CosinePotential()
background = KerrBackground(spin)
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=200.0, step=0.4)

primary, secondary = BoundMode(2, 1, 1), BoundMode(3, 2, 2)
problem_for(modes, potential) = RadialProblem(background, potential, mass_coupling, ModeSet(modes); mesh)

solve_pair(potential) = amplitudes_for_energy(problem_for([primary, secondary], potential), [target_energy, target_energy])
solve_single(mode, potential) = amplitudes_for_energy(problem_for([mode], potential), [target_energy])

@info "linear"
linear = solve_pair(FreePotential())
@info "single-cloud nonlinear"
single_primary = solve_single(primary, potential)
single_secondary = solve_single(secondary, potential)
@info "joint nonlinear"
joint = solve_pair(potential)

profile_of(solution, mode_index) = abs.(profiles(solution)[:, fundamental_channels(solution.problem.mode_set)[mode_index]])

figure = Figure(size=(1400, 650))
axis = Axis(figure[1, 1]; xlabel="r*/M", ylabel="|R_{l,q}|",
            title=@sprintf("Joint-deformation comparison: α = %.2f, a/M = %.2f, E₁ = E₂ ≈ %.1f f_a²M, %d mesh points",
                           mass_coupling, spin, target_energy, length(mesh.rstar)))
colors = (primary=Makie.wong_colors()[1], secondary=Makie.wong_colors()[2])
lines!(axis, mesh.rstar, profile_of(linear, 1); color=colors.primary, linestyle=:dot, label="Primary: linear")
lines!(axis, mesh.rstar, profile_of(single_primary, 1); color=colors.primary, linestyle=:dash, label="Primary: single cloud")
lines!(axis, mesh.rstar, profile_of(joint, 1); color=colors.primary, label="Primary: joint")
lines!(axis, mesh.rstar, profile_of(linear, 2); color=colors.secondary, linestyle=:dot, label="Secondary: linear")
lines!(axis, mesh.rstar, profile_of(single_secondary, 1); color=colors.secondary, linestyle=:dash, label="Secondary: single cloud")
lines!(axis, mesh.rstar, profile_of(joint, 2); color=colors.secondary, label="Secondary: joint")
axislegend(axis; position=:rt, nbanks=2)
save(joinpath(@__DIR__, "profile_comparison.png"), figure)

for (name, solution) in (("linear", linear), ("single primary", single_primary), ("single secondary", single_secondary), ("joint", joint))
    @printf("%-17s amplitudes = %s   frequencies = %s\n", name, solution.amplitudes, mode_frequencies(solution))
end
