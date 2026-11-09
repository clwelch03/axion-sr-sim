# Wave actions N_i = ∂⟨L⟩/∂ω_i of the independent phases (a = 0.9, α = 0.3, cosine, two modes):
#   1. identities: J = Σ m_i N_i, L_E = Σ ω_i L_{N_i}, L_J = Σ m_i L_{N_i}, ∂_A E = Σ ω_i ∂_A N_i, at a linear and a
#      nonlinear point;
#   2. energy versus action balance on a fixed background against A' = γA (linear regime), with the energy flux as
#      coded and with the omitted −Im ω̃ Re(R*Q) term, and the dependence of E − Σ ω N on the outer endpoint;
#   3. the action-based metric work against the Boyer–Lindquist expression and the linear prediction Σ N_i dω_i/dτ,
#      as a function of the inner endpoint.
# Run from the package root:  julia --project=. validation/action_balance.jl
using AxionCloud, Printf, LinearAlgebra
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
azimuthal = [1, 2]

"""Energy flux including the growth term: 4π(r²+a²) Σ_q [ω_q Im(R*Q) − Im ω̃_q Re(R*Q)]."""
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

function identities(solution, label)
    omega = real.(mode_frequencies(solution))
    charges, actions = cloud_charges(solution), phase_actions(solution)
    losses, action_loss = boundary_losses(solution), action_losses(solution)
    directions = amplitude_response(solution)
    charge_jac, action_jac = charge_jacobian(solution; directions), action_jacobian(solution; directions)
    println("== ", label, ": A = ", solution.amplitudes, ", N = ", actions)
    @printf("  (J − Σ m N)/J          = %+.1e\n", (charges.angular_momentum - dot(azimuthal, actions)) / charges.angular_momentum)
    @printf("  (L_E − Σ ω L_N)/L_E    = %+.1e\n", (losses.energy - dot(omega, action_loss)) / losses.energy)
    @printf("  (L_J − Σ m L_N)/L_J    = %+.1e\n", (losses.angular_momentum - dot(azimuthal, action_loss)) / losses.angular_momentum)
    @printf("  (E − Σ ω N)/E          = %+.1e  (outer surface term plus potential energy)\n", (charges.energy - dot(omega, actions)) / charges.energy)
    @printf("  (∂_A E − Σ ω ∂_A N)/∂_A E = (%+.1e, %+.1e)\n", ((charge_jac[1, :] .- vec(transpose(omega) * action_jac)) ./ charge_jac[1, :])...)
    complete_loss = complete_energy_flux(solution, solution.layout.n_nodes) - complete_energy_flux(solution, 1)
    rates = (("energy balance, coded flux", -(charge_jac \ collect(losses))),
             ("energy balance, complete flux", -(charge_jac \ [complete_loss, losses.angular_momentum])),
             ("action balance", -(action_jac \ action_loss)))
    for (name, rate) in rates
        @printf("  %-30s A' = (%.8e, %.8e)   A'/(γA) − 1 = (%+.1e, %+.1e)\n", name, rate...,
                (rate ./ (growth_rates(solution) .* solution.amplitudes) .- 1)...)
    end
    flush(stdout)
end

println("1–2. identities and balances (default mesh r* ∈ [−80, 180])")
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
linear = solve_radial(problem, [1e-6, 1e-7])
identities(linear, "linear")
let solution = linear
    for fraction in exp10.(range(-4, 0; length=9)[2:end])
        solution = solve_radial(problem, [1e-2, 1e-4] .* fraction; guess=solution)
    end
    identities(solution, "nonlinear (A'/(γA) is not expected to be 1 here)")
end

println("\nE − Σ ω N against the outer endpoint (free field, linear point):")
for rstar_max in (180.0, 250.0, 320.0)
    mesh = RadialMesh(background; rstar_min=-80.0, rstar_max, step=0.4)
    solution = solve_radial(RadialProblem(background, FreePotential(), 0.3, mode_set; mesh), [1e-6, 1e-7])
    @printf("  r*_max = %5.0f: (E − Σ ω N)/E = %+.2e\n", rstar_max,
            (cloud_charges(solution).energy - dot(real.(mode_frequencies(solution)), phase_actions(solution))) / cloud_charges(solution).energy)
end

println("\n3. metric work against the inner endpoint (linear point, f_a/M_Pl = 1, b' and j' from the horizon fluxes)")
state = BackgroundState(1.0, background.a)
for rstar_min in (-45.0, -60.0, -80.0, -100.0, -120.0)
    mesh = RadialMesh(background; rstar_min, rstar_max=180.0, step=0.4)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    cache = EvolutionCache(problem; decay_constant_ratio=1.0, last_solution=solution)
    terms = balance_terms(cache, solution, state)
    rate, derivatives = terms.state_rate, terms.derivatives
    actions = AC.reference_actions(phase_actions(solution), state)
    linear_prediction = sum(actions .* (derivatives.frequencies_wrt_mass .* rate.mass .+
                                        derivatives.frequencies_wrt_angular_momentum .* rate.angular_momentum))
    old_work = boyer_lindquist_metric_work(solution, state, rate)
    energy_rate = cloud_charges(solution).energy * rate.mass
    @printf("  r*_min = %6.1f: w/(E b') action-based = %+.6e, linear prediction Σ N dω/dτ = %+.6e, Boyer–Lindquist = %+.3e\n",
            rstar_min, terms.metric_work / energy_rate, linear_prediction / energy_rate, old_work / energy_rate)
    rates_energy = amplitude_rates(terms; balance=:energy)
    rates_energy_work = amplitude_rates(terms; balance=:energy, include_metric_work=true)
    rates_action = amplitude_rates(terms; balance=:action)
    @printf("                  A' energy = %s, energy + w = %s, action = %s\n", rates_energy, rates_energy_work, rates_action)
    flush(stdout)
end
