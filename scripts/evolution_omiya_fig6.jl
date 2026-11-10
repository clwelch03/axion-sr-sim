# Joint (mutually deforming) two-mode evolution at the setup of Omiya et al. 2211.01949, Fig. 6 (their
# non-perturbative primary/secondary calculation overlaid with their perturbative one): a = 0.99, |211⟩ + |322⟩,
# E₁ = 0.27 f_a²M (α = 0.42) or 12 f_a²M (α = 0.3), E₂ = 1e-8 f_a²M, fixed background, followed to ω₁,I t = 10
# (ω₁,I: linear growth rate of |211⟩). Joint solution on a domain that contains the cloud, integrated outer data,
# action balance. Potential: `cosine` (their Eq. 1) or `chiral` (pion-integrated QCD potential, z = m_u/m_d = 2.2/4.79).
# The initial energies are the fundamental-component energies (`amplitudes_for_energy`); their E_i is each cloud's
# energy from its own profile (Eq. 64), which agrees at these amplitudes.
# Run from the package root:  julia --project=. scripts/evolution_omiya_fig6.jl <α> <cosine|chiral>   (submit to Slurm)
using AxionCloud, Printf

mass_coupling = parse(Float64, ARGS[1])
potential_name = ARGS[2]
potential = potential_name == "cosine" ? CosinePotential() :
            potential_name == "chiral" ? ChiralPotential(; quark_mass_ratio=2.2 / 4.79) :
            error("potential must be cosine or chiral, got $potential_name")
initial_energies = [mass_coupling ≈ 0.42 ? 0.27 : mass_coupling ≈ 0.3 ? 12.0 : error("α must be 0.3 or 0.42"), 1e-8]

background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, potential, mass_coupling, mode_set; mesh)
prefix = joinpath(@__DIR__, "..", "results", @sprintf("omiya_fig6_a0.99_alpha%.2f_%s", mass_coupling, potential_name))
@printf("α = %.2f, %s potential, mesh: %d nodes, r* ∈ [%.2f, %.2f]; balance = :action\n", mass_coupling,
        potential_name, length(mesh.rstar), mesh.rstar[1], mesh.rstar[end])

linear_rate = growth_rates(solve_radial(with_potential(problem, FreePotential()), [1e-6, 1e-6]))[1]
seeded = amplitudes_for_energy(problem, initial_energies)
amplitudes = seeded.amplitudes
@printf("seed: A = (%.6e, %.6e), fundamental E = (%.6e, %.6e), ω N = (%.6e, %.6e), E_total = %.6e; ω₁,I = %.6e\n",
        amplitudes..., fundamental_energies(seeded)..., (real.(mode_frequencies(seeded)) .* phase_actions(seeded))...,
        cloud_charges(seeded).energy, linear_rate)
flush(stdout)

span = (0.0, 10.0 / linear_rate)
t = @elapsed result = evolve_fixed_background(problem, amplitudes, span; saveat=span[2] / 200, verbose=true)
println("done: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s"); flush(stdout)

open(prefix * "_diagnostics.csv", "w") do io
    println(io, "# linear_growth_rate_211 = $(linear_rate)")
    println(io, "time,A1,A2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,E,J,L_E,L_J,N1,N2,L_N1,L_N2")
    for d in result.diagnostics
        row = (d.time, d.amplitudes..., real.(d.mode_frequencies)..., imag.(d.mode_frequencies)..., d.charges.energy,
               d.charges.angular_momentum, d.losses.energy, d.losses.angular_momentum, d.actions..., d.action_losses...)
        println(io, join((@sprintf("%.10e", value) for value in row), ","))
    end
end
println("wrote $(prefix)_diagnostics.csv"); flush(stdout)
