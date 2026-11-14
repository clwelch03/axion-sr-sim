# Check of the GW term in balance_terms at one state (needs the gw_power_coefficients keyword of EvolutionCache).
# Run from the package root: julia --project=. validation/gw_rate_check.jl patched
# With the argument "original" only the rates without GW are printed (for comparing before and after the change).
using AxionCloud, Printf
patched = ARGS[1] == "patched"
decay_constant_ratio = 0.0466873
background = KerrBackground(0.9)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
state = BackgroundState(0.99, 0.875 * 0.99^2)   # a hole that has already lost some mass and spin
for (name, potential, target_actions) in (("free", FreePotential(), [28.0, 1e-3]), ("cosine", CosinePotential(), [28.0, 1e-3]))
    problem = RadialProblem(background, potential, mass_coupling, mode_set; mesh)
    current = with_background(problem, KerrBackground(state.angular_momentum / state.mass^2), mass_coupling * state.mass)
    probe = solve_radial(current, [1e-6, 1e-6])
    amplitudes = probe.amplitudes .* sqrt.(target_actions ./ (state.mass^2 .* phase_actions(probe)))
    solution = probe
    for factor in exp.(range(log(1e-3), 0.0; length=20))
        solution = solve_radial(current, amplitudes .* factor; guess=solution)
        solution.converged || error("continuation failed")
    end
    plain_cache = EvolutionCache(problem; decay_constant_ratio, last_solution=solution)
    plain = balance_terms(plain_cache, solution, state)
    plain_rates = amplitude_rates(plain, plain_cache)
    @printf("%s: A = (%.6e, %.6e); without GW: A'/A = (%.16e, %.16e), L_N = (%.16e, %.16e), L_E = %.16e\n", name,
            solution.amplitudes..., (plain_rates ./ solution.amplitudes)..., plain.action_losses..., plain.losses.energy)
    patched || continue
    include(joinpath(@__DIR__, "..", "scripts", "gw_power.jl"))
    table = read_gw_power_table()
    gw_power_coefficients(alpha) = [gw_power_coefficient(table, mode.m, alpha) for mode in mode_set.modes]
    gw_cache = EvolutionCache(problem; decay_constant_ratio, last_solution=solution, gw_power_coefficients)
    gw = balance_terms(gw_cache, solution, state)
    gw_rates = amplitude_rates(gw, gw_cache)
    coefficients = gw_power_coefficients(mass_coupling * state.mass)
    frequencies = real.(mode_frequencies(solution)) ./ state.mass
    actions = state.mass^2 .* phase_actions(solution)
    # free field: N_i ∝ A_i² and the modes are uncoupled, so the GW term shifts A_i'/A_i by −c_i κ ω_i N_i / (2 b²)
    expected_shift = -coefficients .* decay_constant_ratio^2 .* frequencies .* actions ./ (2state.mass^2)
    shift = (gw_rates .- plain_rates) ./ solution.amplitudes
    @printf("  GW shift of A'/A: (%.10e, %.10e); free-field formula: (%.10e, %.10e); ratio − 1: (%.1e, %.1e)\n",
            shift..., expected_shift..., (shift ./ expected_shift .- 1)...)
    cloud_fraction = decay_constant_ratio^2 * frequencies[1] * actions[1] / state.mass
    @printf("  |211⟩: M_cloud/M = %.5f; radiated power κ ΔL_E = %.6e; c (M_cloud/M)² (211 + 322) = %.6e\n", cloud_fraction,
            decay_constant_ratio^2 * (gw.losses.energy - plain.losses.energy),
            sum(coefficients .* (decay_constant_ratio^2 .* frequencies .* actions ./ state.mass) .^ 2))
    @printf("  ΔL_J − Σ m_i ΔL_N,i = %.1e; hole rates unchanged: %s\n",
            gw.losses.angular_momentum - plain.losses.angular_momentum - sum([1, 2] .* gw.gw_action_losses),
            gw.state_rate == plain.state_rate)
    energy_cache = EvolutionCache(problem; decay_constant_ratio, last_solution=solution, gw_power_coefficients, balance=:energy)
    plain_energy_cache = EvolutionCache(problem; decay_constant_ratio, last_solution=solution, balance=:energy)
    energy_shift = (amplitude_rates(gw, energy_cache) .- amplitude_rates(plain, plain_energy_cache)) ./ solution.amplitudes
    @printf("  energy balance: GW shift of A'/A = (%.10e, %.10e); relative to action balance − 1: (%.1e, %.1e)\n",
            energy_shift..., (energy_shift ./ shift .- 1)...)
end
