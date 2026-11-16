# Issues G and H (reviewer): what the wave actions contain and how they depend on the outer endpoint, and whether
# ∂_A E = Σ ω_i ∂_A N_i holds, as a function of the outer and inner endpoints (a = 0.9, α = 0.3, cosine potential,
# fixed background). Usage:
#   julia --project=. validation/domain_actions_check.jl [outer|inner] [amplitude1 amplitude2]
# "outer" varies r*_max at inner endpoint −60 (plus the default estimated mesh); "inner" varies r*_min at r*_max = 500.
using AxionCloud, Printf, LinearAlgebra
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
alpha = 0.3
amplitude_target = length(ARGS) >= 3 ? parse.(Float64, ARGS[2:3]) : [1e-2, 1e-3]
which = length(ARGS) >= 1 ? ARGS[1] : "outer"
radiating = findfirst(c -> c.q == [-1, 2], mode_set.components)
tail_cuts = (200.0, 300.0)   # tortoise radii beyond which the (-1,2) charge is removed in the "tail removed" variants

function complete_energy_flux(solution, node)
    problem, layout = solution.problem, solution.layout
    state = AC.physical_state(layout, solution.unknowns, solution.scale, node)
    frequencies = AC.component_frequencies(problem.mode_set, mode_frequencies(solution))
    total = 0.0
    for (index, channels) in enumerate(problem.mode_set.component_channels), channel in channels
        product = conj(state.profiles[channel]) * state.slopes[channel]
        total += real(frequencies[index]) * imag(product) - imag(frequencies[index]) * real(product)
    end
    return 4pi * AC.r2a2(problem.background, problem.mesh.offsets[node]) * total
end

"""Charges n_q of every component over the full mesh, then n_radiating over r* ≤ each cut (so n_tail = full − cut)."""
function component_charge_observable(problem, layout, unknowns, scale)
    mesh = problem.mesh
    frequencies = AC.real_component_frequencies(problem, layout, unknowns)
    weights = AC.simpson_weights(mesh.rstar)
    totals = zero(AC.component_action_densities(problem, mesh.offsets[1], AC.physical_state(layout, unknowns, scale, 1).profiles, frequencies))
    inside = zeros(eltype(unknowns), length(tail_cuts))
    for node in 1:layout.n_nodes
        density = AC.component_action_densities(problem, mesh.offsets[node],
                      AC.physical_state(layout, unknowns, scale, node).profiles, frequencies)
        totals = totals .+ weights[node] .* density
        for (cut_index, cut) in enumerate(tail_cuts)
            mesh.rstar[node] <= cut && (inside[cut_index] += weights[node] * density[radiating])
        end
    end
    return vcat(totals, inside)   # (weights are not truncated at the cut; error ≤ one Simpson cell)
end

function make_mesh(label)
    label == :default && return RadialMesh(background, alpha, mode_set)
    rstar_min, rstar_max = label
    return RadialMesh(background; rstar_min, rstar_max, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step=3.2)
end

function solve_nonlinear(mesh, target)
    problem = RadialProblem(background, CosinePotential(), alpha, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    maximum(target) < 1e-4 && return solve_radial(problem, target; guess=solution)
    # AXION_CONTINUATION_STEPS (default 7): number of log-spaced amplitude steps from 1e-4 of the target to the target
    for fraction in exp10.(range(-4, 0; length=parse(Int, get(ENV, "AXION_CONTINUATION_STEPS", "7")))[2:end])
        solution = solve_radial(problem, target .* fraction; guess=solution)
        solution.converged || @warn "non-converged continuation step" fraction solution.residual_norm
    end
    return solution
end

function report(label, solution)
    problem = solution.problem
    mesh = problem.mesh
    omega = real.(mode_frequencies(solution))
    weights_matrix = AC.phase_weights(mode_set)
    directions = amplitude_response(solution)
    energy_jac = charge_jacobian(solution; directions)
    action_jac = action_jacobian(solution; directions)
    comp_jac = AC.observable_jacobian(solution, component_charge_observable; directions)
    comp_values = component_charge_observable(problem, solution.layout, solution.unknowns, solution.scale)
    n_components = length(mode_set.components)
    actions, action_loss = phase_actions(solution), action_losses(solution)
    losses = boundary_losses(solution)
    charges = cloud_charges(solution)
    rate_action = -(action_jac \ action_loss)
    rate_energy = -(energy_jac \ collect(losses))
    complete_loss = complete_energy_flux(solution, solution.layout.n_nodes) - complete_energy_flux(solution, 1)
    rate_energy_complete = -(energy_jac \ [complete_loss, losses.angular_momentum])
    D = vec(energy_jac[1, :]) .- vec(transpose(omega) * action_jac)
    dE_dtau = dot(energy_jac[1, :], rate_action)
    wc = weights_matrix[:, radiating]
    tail_full(cut_index) = comp_values[radiating] - comp_values[n_components + cut_index]
    println("---- ", label, "  r* ∈ [", @sprintf("%.1f, %.1f", mesh.rstar[1], mesh.rstar[end]), "], nodes ", length(mesh.rstar),
            ", newton iterations ", solution.iterations, ", converged ", solution.converged)
    @printf("  A = %s  ω = %s\n", solution.amplitudes, mode_frequencies(solution))
    @printf("  N = (%.8e, %.8e)   E = %.8e   (E − ΣωN)/E = %+.2e\n", actions..., charges.energy,
            (charges.energy - dot(omega, actions)) / charges.energy)
    @printf("  n_q:"); for (i, c) in enumerate(mode_set.components); @printf("  %s: %.6e", string(c.q), comp_values[i]); end; println()
    @printf("  component_charges energy of (-1,2) = %.6e\n", component_charges(solution, radiating).energy)
    @printf("  tail of n_(-1,2) beyond r*=%.0f: %.4e (%.2e of N_2);  beyond %.0f: %.4e\n", tail_cuts[1], tail_full(1),
            2tail_full(1) / actions[2], tail_cuts[2], tail_full(2))
    @printf("  dN/dA = [[%.6e %.6e]; [%.6e %.6e]]\n", transpose(action_jac)[1, :]..., transpose(action_jac)[2, :]...)
    @printf("  d n_(-1,2)/dA = (%.6e, %.6e)\n", comp_jac[radiating, :]...)
    @printf("  L_N = (%.6e, %.6e)   L_E = %.6e  L_J = %.6e  complete L_E = %.6e\n", action_loss..., losses.energy, losses.angular_momentum, complete_loss)
    @printf("  Adot action          = (%.8e, %.8e)\n", rate_action...)
    @printf("  Adot energy (coded)  = (%.8e, %.8e)  ratio−1 vs action (%+.1e, %+.1e)\n", rate_energy..., (rate_energy ./ rate_action .- 1)...)
    @printf("  Adot energy (complete flux) ratio−1 vs action (%+.1e, %+.1e)\n", (rate_energy_complete ./ rate_action .- 1)...)
    @printf("  dE/dA = (%.6e, %.6e)   D = dE/dA − Σω dN/dA = (%+.3e, %+.3e)   D/|dE/dA| = (%+.2e, %+.2e)\n",
            energy_jac[1, :]..., D..., (D ./ abs.(energy_jac[1, :]))...)
    @printf("  Σ D Adot / (dE/dτ) = %+.3e   (dE/dτ = %.6e)\n", dot(D, rate_action) / dE_dtau, dE_dtau)
    # tail-removed variants: N_i → N_i − W[i,c] (n_c − n_c(r* ≤ cut)), Jacobian likewise; losses unchanged
    for (name, removed_value, removed_jac) in (
            ("tail beyond r*=200 removed", tail_full(1), comp_jac[radiating, :] .- comp_jac[n_components + 1, :]),
            ("tail beyond r*=300 removed", tail_full(2), comp_jac[radiating, :] .- comp_jac[n_components + 2, :]),
            ("(-1,2) removed entirely", comp_values[radiating], comp_jac[radiating, :]))
        jac_variant = action_jac .- wc * transpose(removed_jac)
        rate_variant = -(jac_variant \ action_loss)
        @printf("  Adot, %-28s = (%.8e, %.8e)  rel. change (%+.2e, %+.2e)\n", name, rate_variant...,
                (rate_variant ./ rate_action .- 1)...)
    end
    # tail shape: density of n_(-1,2) per unit r* at a few radii, and k
    frequencies = AC.component_frequencies(mode_set, mode_frequencies(solution))
    wave_number = outer_wavenumber(frequencies[radiating], alpha)
    omega_real = AC.real_component_frequencies(problem, solution.layout, solution.unknowns)
    @printf("  k(-1,2) = %.5e %+.3ei;  dn/dr* of (-1,2):", real(wave_number), imag(wave_number))
    for rstar_probe in (100.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0)
        rstar_probe <= mesh.rstar[end] || continue
        node = argmin(abs.(mesh.rstar .- rstar_probe))
        density = AC.component_action_densities(problem, mesh.offsets[node],
                      AC.physical_state(solution.layout, solution.unknowns, solution.scale, node).profiles, omega_real)[radiating]
        @printf("  r*=%.0f: %.3e", mesh.rstar[node], density)
    end
    println(); flush(stdout)
end

# optional 4th argument: run only the i-th mesh of the list (so the meshes can run as parallel processes)
pick(list) = length(ARGS) >= 4 ? [list[parse(Int, ARGS[4])]] : list
if which == "outer"
    # AXION_OUTER_ENDPOINT, if set, replaces the list by that single outer tortoise endpoint
    outer_labels = haskey(ENV, "AXION_OUTER_ENDPOINT") ? [(-60.0, parse(Float64, ENV["AXION_OUTER_ENDPOINT"]))] :
                   pick([(-60.0, 500.0), (-60.0, 700.0), (-60.0, 1000.0), :default])
    for label in outer_labels
        t = @elapsed solution = solve_nonlinear(make_mesh(label), amplitude_target)
        @printf("[%s solve time %.0f s]\n", string(label), t)
        report(string(label), solution)
    end
else
    for rstar_min in pick([-40.0, -60.0, -80.0])
        t = @elapsed solution = solve_nonlinear(make_mesh((rstar_min, 500.0)), amplitude_target)
        @printf("[inner %.0f solve time %.0f s]\n", rstar_min, t)
        report("inner $(rstar_min)", solution)
    end
end
