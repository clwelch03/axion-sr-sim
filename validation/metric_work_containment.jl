# Action-based metric work against its linear-regime prediction Σ N_i dω_i/dτ as the outer endpoint moves out
# (a = 0.9, α = 0.3, linear point, f_a/M_Pl = 1, inner endpoint −60). The two must agree once the cloud is contained,
# because they differ exactly by ∂_b,j of the surface term E − Σ ω N.
# Run from the package root:  julia --project=. validation/metric_work_containment.jl
using AxionCloud, Printf, LinearAlgebra
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
state = BackgroundState(1.0, background.a)
for rstar_max in (180.0, 250.0, 320.0, 400.0)
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max, step=0.4)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    terms = balance_terms(EvolutionCache(problem; decay_constant_ratio=1.0, last_solution=solution), solution, state)
    rate, derivatives = terms.state_rate, terms.derivatives
    actions = AC.reference_actions(phase_actions(solution), state)
    prediction = sum(actions .* (derivatives.frequencies_wrt_mass .* rate.mass .+ derivatives.frequencies_wrt_angular_momentum .* rate.angular_momentum))
    energy_rate = cloud_charges(solution).energy * rate.mass
    surface = (cloud_charges(solution).energy - dot(real.(mode_frequencies(solution)), phase_actions(solution))) / cloud_charges(solution).energy
    rates_energy, rates_work, rates_action = amplitude_rates(terms), amplitude_rates(terms; include_metric_work=true), amplitude_rates(terms; balance=:action)
    @printf("r*_max = %5.0f: (E − ΣωN)/E = %+.1e  w/(E b') = %+.6e  Σ N dω/dτ/(E b') = %+.6e  A'_action/A'_energy − 1 = (%+.1e, %+.1e)  with w: (%+.1e, %+.1e)\n",
            rstar_max, surface, terms.metric_work / energy_rate, prediction / energy_rate,
            (rates_action ./ rates_energy .- 1)..., (rates_action ./ rates_work .- 1)...)
    flush(stdout)
end
