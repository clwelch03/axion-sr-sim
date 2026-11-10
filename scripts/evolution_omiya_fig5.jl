# Joint (mutually deforming) two-mode evolution at the parameters of Omiya et al. 2211.01949, Fig. 5:
# a = 0.99, cosine potential, |211⟩ and |322⟩ seeded with E₁ = E₂ = 1e-70 f_a² M, fixed background, followed to
# ω₁,I t = 87 (ω₁,I: linear growth rate of |211⟩). Their figure uses their primary/secondary approximation; this run
# uses the joint solution on a domain that contains the cloud, integrated outer data and the action balance.
# Run from the package root:  julia --project=. scripts/evolution_omiya_fig5.jl <α>   (α = 0.3 or 0.42; submit to Slurm)
using AxionCloud, Printf

mass_coupling = parse(Float64, ARGS[1])
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
prefix = joinpath(@__DIR__, "..", "results", @sprintf("omiya_fig5_a0.99_alpha%.2f", mass_coupling))
@printf("α = %.2f, mesh: %d nodes, r* ∈ [%.2f, %.2f]; balance = :action\n", mass_coupling, length(mesh.rstar),
        mesh.rstar[1], mesh.rstar[end])

# Seed: E_i = 1e-70 each. In the linear regime E_i = ω_i N_i and N_i ∝ A_i², so one probe solve fixes the amplitudes.
initial_energy = 1e-70
probe = solve_radial(problem, [1e-6, 1e-6])
target_actions = initial_energy ./ real.(mode_frequencies(probe))
amplitudes = probe.amplitudes .* sqrt.(target_actions ./ phase_actions(probe))
seeded = solve_radial(problem, amplitudes; guess=probe)
growth_rate = growth_rates(seeded)[1]
@printf("seed: A = (%.6e, %.6e), N = (%.6e, %.6e), ω N = (%.6e, %.6e); ω₁,I = %.6e\n", amplitudes...,
        phase_actions(seeded)..., (real.(mode_frequencies(seeded)) .* phase_actions(seeded))..., growth_rate)
flush(stdout)

span = (0.0, 87.0 / growth_rate)
t = @elapsed result = evolve_fixed_background(problem, amplitudes, span; saveat=span[2] / 200, verbose=true)
println("done: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s"); flush(stdout)

open(prefix * "_diagnostics.csv", "w") do io
    println(io, "# growth_rate_211 = $(growth_rate)")
    println(io, "time,A1,A2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,E,J,L_E,L_J,N1,N2,L_N1,L_N2")
    for d in result.diagnostics
        row = (d.time, d.amplitudes..., real.(d.mode_frequencies)..., imag.(d.mode_frequencies)..., d.charges.energy,
               d.charges.angular_momentum, d.losses.energy, d.losses.angular_momentum, d.actions..., d.action_losses...)
        println(io, join((@sprintf("%.10e", value) for value in row), ","))
    end
end
println("wrote $(prefix)_diagnostics.csv"); flush(stdout)
