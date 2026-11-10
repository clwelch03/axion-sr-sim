# Checkpoint for the Omiya-versus-joint rate map: at fixed |211⟩ amplitude A₁, how do the per-component fluxes of the
# joint solution depend on the |322⟩ amplitude A₂?
#
# Omiya et al. 2211.01949 Sec. III (Eqs. 47–52) keep |211⟩ fully nonlinear but independent of A₂, take |322⟩ linear on
# that background (φ₂ = A₂ φ̂₂(A₁)), and source the radiative modes from φ₂ to the order it first appears. So in their
# model the flux carried by a Fourier component q = (q₁, q₂) (frequency q₁ω₁ + q₂ω₂) is the leading power of A₂ in that
# component, c_q(A₁) |A₂|^(2|q₂|): A₂-independent for (1,0), ∝ |A₂|² for (0,1) and the m = 0 wave (2,−1), ∝ |A₂|⁴ for
# the m = 3 wave (−1,2). Here c_q is read off a joint solve at small A₂ and extrapolated; the joint solve at the full A₂
# is our rate. Part 1 checks that the small-A₂ coefficients are clean (constant as A₂ → 0); part 2 shows the departure.
#
# Run from the package root:  julia --project=. scripts/omiya_flux_extrapolation.jl [α] [E₂₁₁]   (defaults 0.3, 300)
using AxionCloud, Printf

mass_coupling = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 0.3
primary_energy = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 300.0

background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
components = mode_set.components
leading_power = [2abs(component.q[2]) for component in components]

"""Outward action flux of each Fourier component (unweighted by phase) at mesh `node`; same lines as `action_flux`."""
function component_action_fluxes(solution, node)
    state = AxionCloud.physical_state(solution.layout, solution.unknowns, solution.scale, node)
    wronskian = imag.(conj.(state.profiles) .* state.slopes)
    prefactor = 4pi * AxionCloud.r2a2(problem.background, problem.mesh.offsets[node])
    [prefactor * sum(wronskian[channels]) for channels in mode_set.component_channels]
end
"""Net loss of each component, outer minus inner flux (so that summing with phase weights gives `action_losses`)."""
component_losses(solution) = component_action_fluxes(solution, solution.layout.n_nodes) -
                             component_action_fluxes(solution, 1)

# Fix A₁ from a calibration at E₂₁₁ = primary_energy with a tiny secondary, then scan A₂ at that A₁.
calibrated = amplitudes_for_energy(problem, [primary_energy, 1e-6])
primary_amplitude, reference_secondary = calibrated.amplitudes
@printf("α = %.2f, a = 0.99, cosine; A₁ = %.6e fixed (E₂₁₁ = %.4g at E₃₂₂ = 1e-6)\n", mass_coupling,
        primary_amplitude, primary_energy)
println("components q: ", join(("(" * join(component.q, ",") * ")" for component in components), "  "))
flush(stdout)

secondary_energies = [1e-6, 1e-4, 1e-2, 1.0, 10.0, 30.0, 100.0, 300.0]
rows = []
guess = calibrated
for secondary_energy in secondary_energies
    secondary_amplitude = reference_secondary * sqrt(secondary_energy / 1e-6)
    elapsed = @elapsed solution = solve_radial(problem, [primary_amplitude, secondary_amplitude]; guess)
    solution.converged || (@printf("E₃₂₂ = %g: not converged, stopping\n", secondary_energy); break)
    global guess = solution
    push!(rows, (; secondary_amplitude, solution, energies=fundamental_energies(solution), frequencies=mode_frequencies(solution),
                 losses=component_losses(solution), action_losses=action_losses(solution)))
    @printf("solved E₃₂₂ ≈ %-7g  A₂ = %.4e  fundamental E = (%.4g, %.4g)  Re ω = (%.8f, %.8f)  [%.0f s]\n",
            secondary_energy, secondary_amplitude, rows[end].energies..., real.(rows[end].frequencies)..., elapsed)
    flush(stdout)
end

# Part 1: coefficient c_q = loss_q / A₂^(2|q₂|), relative to the reference solve. Constant ⇒ clean extrapolation.
# Reference: E₃₂₂ = 1e-2. At 1e-6 the |A₂|⁴ and |A₂|⁶ components are below the solver's resolution (first run of this
# script: the (−1,2) coefficient changed by a factor 20 between E₃₂₂ = 1e-6 and 1e-4, then stayed fixed to 1e-5).
reference_index = findfirst(==(1e-2), secondary_energies)
reference = rows[reference_index]
println("\nPart 1: loss_q / A₂^p relative to the E₃₂₂ = 1e-2 solve (p = 2|q₂|), and the reference value itself")
println("  raw component losses at the reference: ", join((@sprintf("%+.3e", loss) for loss in reference.losses), " "))
coefficients = [reference.losses[k] / reference.secondary_amplitude^leading_power[k] for k in eachindex(components)]
for (k, component) in enumerate(components)
    @printf("  q = (%2d,%2d) p = %d  c_q = %+.6e   ratios:", component.q..., leading_power[k], coefficients[k])
    for row in rows
        @printf(" %+.6f", row.losses[k] / row.secondary_amplitude^leading_power[k] / coefficients[k])
    end
    println()
end

# Part 2: joint versus extrapolated (Omiya-type) phase action losses L_N, and per component.
weights = AxionCloud.phase_weights(mode_set)
println("\nPart 2: action losses L_N (dN_i/dt = −L_N), joint vs leading-power extrapolation")
for row in rows
    extrapolated = coefficients .* row.secondary_amplitude .^ leading_power
    omiya_losses = weights * extrapolated
    @printf("  E₃₂₂ = %-9.3g joint L_N = (%+.5e, %+.5e)  extrapolated = (%+.5e, %+.5e)  joint/extrap − 1 = (%+.3e, %+.3e)\n",
            row.energies[2], row.action_losses..., omiya_losses..., (row.action_losses ./ omiya_losses .- 1)...)
end
println("\nper-component joint/extrapolated − 1 (columns: E₃₂₂ as above)")
for (k, component) in enumerate(components)
    @printf("  q = (%2d,%2d):", component.q...)
    for row in rows
        @printf(" %+.3e", row.losses[k] / (coefficients[k] * row.secondary_amplitude^leading_power[k]) - 1)
    end
    println()
end
flush(stdout)

# Part 3: the same comparison at equal occupations (N₂₁₁, N₃₂₂), which is what the evolution law dN/dt = −L_N(N)
# compares. Omiya-type rates at (N₁, N₂): find the A₁ whose solve with a tiny secondary (E₃₂₂ ≈ 1e-2) has phase action
# N₁, then scale each component loss by (N₂ / N₂,ref)^|q₂| (in their model N₂ ∝ |A₂|² at fixed A₁).
"""Small-secondary solve whose |211⟩ phase action is `primary_action`, by fixed-point iteration A₁ ← A₁ √(N₁/N₁(A₁))."""
function primary_matched_reference(primary_action, start)
    solution = start
    primary, secondary = start.amplitudes
    for _ in 1:20
        primary *= sqrt(primary_action / phase_actions(solution)[1])
        solution = solve_radial(problem, [primary, secondary]; guess=solution)
        abs(phase_actions(solution)[1] / primary_action - 1) < 1e-9 && return solution
    end
    error("A₁ matching did not converge for N₁ = $primary_action")
end

reference_solution = reference.solution
println("\nPart 3: action losses at equal occupations, joint vs Omiya-type (matched N₂₁₁, |322⟩ scaled as N₃₂₂^|q₂|)")
for row_index in eachindex(rows)
    row = rows[row_index]
    actions = phase_actions(row.solution)
    matched = primary_matched_reference(actions[1], reference_solution)
    action_ratio = actions[2] / phase_actions(matched)[2]
    omiya_components = component_losses(matched) .* action_ratio .^ abs.(getindex.(getfield.(components, :q), 2))
    omiya_losses = weights * omiya_components
    @printf("  N = (%.5e, %.5e)  matched A₁ = %.6e (fixed-A₁ value %.6e)\n", actions..., matched.amplitudes[1],
            primary_amplitude)
    @printf("      joint L_N = (%+.5e, %+.5e)  Omiya-type = (%+.5e, %+.5e)  joint/Omiya − 1 = (%+.3e, %+.3e)\n",
            row.action_losses..., omiya_losses..., (row.action_losses ./ omiya_losses .- 1)...)
    @printf("      per component joint/Omiya − 1:")
    for k in eachindex(components)
        @printf(" (%d,%d) %+.3e", components[k].q..., row.losses[k] / omiya_components[k] - 1)
    end
    println(); flush(stdout)
end
