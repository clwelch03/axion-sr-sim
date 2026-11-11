# Joint-deformation comparison of |R_{l,q}| and of the energy per unit tortoise coordinate dE/dr* of the primary
# |211⟩ and secondary |322⟩ at equal per-cloud energy,
# a = 0.99, α = 0.3, with three models:
#   Linear – free field (V = ½μ²φ², no self-interaction), scaled to the target energy; labelled "Free field" in plots;
#   Omiya  – Omiya et al. 2211.01949 Sec. III: primary is the nonlinear single cloud at E₁ (independent of A₂);
#            secondary is linear on that primary, i.e. a solve with the single-cloud A₁ and a tiny A₂ (E₂ ≈ 1e-2),
#            scaled by √(E₂/E₂,ref) (its energy is ∝ |A₂|² at leading order);
#   Joint  – both clouds nonlinear and mutually deforming.
# Energies are the free-field energies of the fundamental components (`fundamental_energies`), in f_a²M; this is our
# energy measure, not Omiya's. All solves share the mesh sized by RadialMesh(background, α, mode set).
#
# The caption states the physical reach: the spin room from a = 0.99 to the |211⟩ threshold Ω_H = ω₂₁₁ bounds the
# cloud angular momentum, Δa = J (f_a/M_Pl)² (J in f_a²M², M_Pl non-reduced), ignoring the drift of α with M.
#
# Run from the package root:
#   julia --project=scripts scripts/profile_comparison_regimes.jl <label> <cosine|chiral> <energy>
# Writes results/profile_comparison_<label>.csv; plot with scripts/plot_profile_comparison_regimes.jl.
using AxionCloud, Printf

label, potential_name, target_energy = ARGS[1], ARGS[2], parse(Float64, ARGS[3])
potential = Dict("cosine" => CosinePotential(), "chiral" => ChiralPotential())[potential_name]
spin, mass_coupling = 0.99, 0.3
reference_secondary_energy = 1e-2

background = KerrBackground(spin)
primary, secondary = BoundMode(2, 1, 1), BoundMode(3, 2, 2)
pair = ModeSet([primary, secondary])
mesh = RadialMesh(background, mass_coupling, pair)
@printf("%s: %s potential, E = %g, mesh %d nodes, r* ∈ [%.1f, %.1f]\n", label, potential_name, target_energy,
        length(mesh.rstar), mesh.rstar[1], mesh.rstar[end])
problem_for(modes, potential) = RadialProblem(background, potential, mass_coupling, ModeSet(modes); mesh)
pair_problem = problem_for([primary, secondary], potential)

timed(name, f) = (t = time(); s = f(); @printf("  %-22s %.0f s, converged = %s, |F| = %.2e\n", name, time() - t,
                                               s.converged, s.residual_norm); flush(stdout); s)
linear = timed("linear pair", () -> amplitudes_for_energy(problem_for([primary, secondary], FreePotential()),
                                                          [target_energy, target_energy]))
single_primary = timed("single primary", () -> amplitudes_for_energy(problem_for([primary], potential), [target_energy]))
# Warm start for the reference: the pair calibrated to (E₁, E₂,ref), whose A₁ differs from the single-cloud one by
# O(E₂,ref); a cold solve at the nonlinear A₁ from the linear guess does not converge at E = 650.
small_secondary_pair = timed("small-secondary pair", () -> amplitudes_for_energy(pair_problem,
                                                                                 [target_energy, reference_secondary_energy];
                                                                                 guess=linear))
omiya_reference = timed("Omiya reference", () -> solve_radial(pair_problem,
                                                              [single_primary.amplitudes[1], small_secondary_pair.amplitudes[2]];
                                                              guess=small_secondary_pair))
omiya_reference.converged || error("Omiya reference solve did not converge")
joint = timed("joint pair", () -> amplitudes_for_energy(pair_problem, [target_energy, target_energy]; guess=linear))
omiya_scale = sqrt(target_energy / fundamental_energies(omiya_reference)[2])

"""
    fundamental_energy_densities(solution)

Matrix `[node, mode]` of dE_i/dr* (in f_a²): the free-field energy per unit tortoise coordinate of each bound
mode's fundamental component alone, i.e. the integrand of `component_charges` before the Simpson weight, so that
integrating a column over r* gives `fundamental_energies(solution)[i]`.
"""
function fundamental_energy_densities(solution)
    problem, layout, scale = solution.problem, solution.layout, solution.scale
    mode_set, mesh = problem.mode_set, problem.mesh
    frequencies = AxionCloud.real_component_frequencies(problem, layout, solution.unknowns)
    fundamental_components = [mode_set.channel_component[channel] for channel in fundamental_channels(mode_set)]
    densities = zeros(layout.n_nodes, length(fundamental_components))
    for node in 1:layout.n_nodes
        state = AxionCloud.physical_state(layout, solution.unknowns, scale, node)
        for (mode_index, component) in enumerate(fundamental_components)
            channels = mode_set.component_channels[component]
            only_this, only_this_slopes = zero(state.profiles), zero(state.slopes)
            only_this[channels], only_this_slopes[channels] = state.profiles[channels], state.slopes[channels]
            bilinears = AxionCloud.field_bilinears(problem, mesh.offsets[node], only_this, only_this_slopes,
                                                   frequencies; potential=FreePotential())
            densities[node, mode_index] = AxionCloud.charge_densities(problem, mesh.offsets[node], bilinears).energy
        end
    end
    # Checkpoint: the densities integrate back to the calibrated energies.
    integrated = vec(sum(AxionCloud.simpson_weights(mesh.rstar) .* densities; dims=1))
    @printf("  dE/dr* integrates to %s (fundamental_energies %s)\n", integrated, fundamental_energies(solution))
    isapprox(integrated, fundamental_energies(solution); rtol=1e-10) || error("energy densities do not integrate to E")
    return densities
end

channels = fundamental_channels(pair)
single_channel = fundamental_channels(ModeSet([primary]))[1]
linear_density, single_density = fundamental_energy_densities(linear), fundamental_energy_densities(single_primary)
omiya_density, joint_density = fundamental_energy_densities(omiya_reference), fundamental_energy_densities(joint)
# Energy is quadratic in the profiles, so the √(E₂/E₂,ref) profile scaling is a factor E₂/E₂,ref on the density.
densities = (energy_density_linear_primary=linear_density[:, 1],
             energy_density_omiya_primary=single_density[:, 1],
             energy_density_joint_primary=joint_density[:, 1],
             energy_density_linear_secondary=linear_density[:, 2],
             energy_density_omiya_secondary=omiya_scale^2 .* omiya_density[:, 2],
             energy_density_joint_secondary=joint_density[:, 2])
curves = merge((linear_primary=abs.(profiles(linear)[:, channels[1]]),
          omiya_primary=abs.(profiles(single_primary)[:, single_channel]),
          joint_primary=abs.(profiles(joint)[:, channels[1]]),
          linear_secondary=abs.(profiles(linear)[:, channels[2]]),
          omiya_secondary=omiya_scale .* abs.(profiles(omiya_reference)[:, channels[2]]),
          joint_secondary=abs.(profiles(joint)[:, channels[2]])), densities)

# Physical reach.
omega_211 = real(mode_frequencies(joint)[1])
threshold_spin = 4omega_211 / (1 + 4omega_211^2)          # Ω_H = a/(2r₊) = ω₂₁₁
spin_room = spin - threshold_spin
joint_charges = cloud_charges(joint)
largest_ratio = sqrt(spin_room / joint_charges.angular_momentum)   # f_a/M_Pl at which the clouds use all the room
qcd_decay_constant_gev(mass_ev) = 5.7e-6 / mass_ev * 1e12        # m_a = 5.7 μeV (1e12 GeV / f_a)
planck_mass_gev = 1.2209e19
solar_masses(mass_ev) = mass_coupling / (7.4803e9 * mass_ev)    # α = G M m_a / (ħ c³)
qcd_mass_ev = 5.7e-6 * 1e12 / (largest_ratio * planck_mass_gev)  # QCD axion with f_a = largest f_a
reach = @sprintf("Spin room a: %.2f → %.3f with J = %.1f fₐ²M² needs fₐ/M_Pl ≤ %.3g (fₐ ≤ %.2g GeV); QCD at that fₐ: mₐ = %.2g eV, M = %.2g M☉",
                 spin, threshold_spin, joint_charges.angular_momentum, largest_ratio, largest_ratio * planck_mass_gev,
                 qcd_mass_ev, solar_masses(qcd_mass_ev))
println(reach)
for (name, solution) in (("linear", linear), ("single primary", single_primary), ("Omiya reference", omiya_reference),
                         ("joint", joint))
    @printf("  %-16s A = %s  E_fund = %s  ω = %s\n", name, solution.amplitudes, fundamental_energies(solution),
            real.(mode_frequencies(solution)))
end

results_path = joinpath(@__DIR__, "..", "results", "profile_comparison_$(label).csv")
open(results_path, "w") do io
    println(io, "# ", reach)
    @printf(io, "# angular_momentum = %.10g, omega_211 = %.10g, threshold_spin = %.10g, spin = %g, alpha = %g, energy = %g, potential = %s\n",
            joint_charges.angular_momentum, omega_211, threshold_spin, spin, mass_coupling, target_energy, potential_name)
    println(io, "r_star_over_M,", join(string.(keys(curves)), ","))
    for node in eachindex(mesh.rstar)
        println(io, mesh.rstar[node], ",", join((curve[node] for curve in curves), ","))
    end
end
println(results_path)
