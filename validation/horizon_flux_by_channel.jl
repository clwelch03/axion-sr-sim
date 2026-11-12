# Why does the hole spin up in the a₀ = 0.9, f_a/M_Pl = 0.047 restart run (job 61110567) once |322⟩ grows?
# Split the horizon (inner-boundary) flux by channel at saved states of that run and check, channel by channel:
#   1. sign: the energy flux into the hole of a channel with frequency ω_c and azimuthal number m_c has the sign of
#      ω_c (ω_c − m_c Ω_H) (absorbed if ω_c > m_c Ω_H, extracted if superradiant). A wrong sign for one channel (e.g.
#      the |322⟩ channel) would not show in energy or angular-momentum conservation, which the hole update satisfies by
#      construction (it uses the same flux).
#   2. J/E = m_c/ω_c per channel (by construction in outward_flux; checked anyway).
#   3. the channel sum equals outward_flux(solution, 1).
# It also gives each channel's contribution to da/dt of the hole, with
#   b' = −κ F_E,   j' = −κ b F_J   (reference units, F = outward flux at the horizon in current units, κ = (f_a/M_Pl)²),
#   a = j/b²  ⇒  da/dt = j'/b² − 2a b'/b.
# States are rows of results/spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart_solution.csv (t, A1, A2, b, j); each
# radial solution is rebuilt by amplitude continuation from the linear regime, as in scripts/evolution_spindown.jl.
# Run from the package root:  julia --project=. validation/horizon_flux_by_channel.jl   (a few minutes)
using AxionCloud, Printf

const C = AxionCloud
decay_constant_ratio = 0.0466873
coupling_to_hole = decay_constant_ratio^2
initial_spin, initial_mass_coupling = 0.9, 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
initial_background = KerrBackground(initial_spin)
template = RadialProblem(initial_background, CosinePotential(), initial_mass_coupling, mode_set;
                         mesh=RadialMesh(initial_background, initial_mass_coupling, mode_set))
problem_at(mass, angular_momentum) =
    with_background(template, KerrBackground(angular_momentum / mass^2), initial_mass_coupling * mass)

solution_csv = joinpath(@__DIR__, "..", "results", "spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart_solution.csv")
rows = [parse.(Float64, split(line, ',')) for line in readlines(solution_csv)[2:end]]
target_times = [1.5e9, 2.1e10, 2.9e10, 3.4e10, 5.1e10, 9.0e10]

for target_time in target_times
    time, amplitude_211, amplitude_322, mass, angular_momentum = rows[argmin([abs(row[1] - target_time) for row in rows])]
    problem = problem_at(mass, angular_momentum)
    amplitudes = [amplitude_211, amplitude_322]
    solution = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=30))
        solution = solve_radial(problem, amplitudes .* factor; guess=solution)
        solution.converged || error("continuation failed at t = $time, factor $factor")
    end
    background = problem.background
    spin = angular_momentum / mass^2
    horizon_angular_velocity = background.horizon_angular_velocity
    layout, unknowns, scale = solution.layout, solution.unknowns, solution.scale
    state = C.physical_state(layout, unknowns, scale, 1)
    component_frequencies = C.real_component_frequencies(problem, layout, unknowns)
    prefactor = 4pi * C.r2a2(background, problem.mesh.offsets[1])
    actions = phase_actions(solution) .* mass^2
    @printf("\nt = %.3e  a = %.6f  M/M₀ = %.6f  𝒩 = (%.4g, %.4g)  Ω_H = %.6f (current units)  Im ω = (%.3e, %.3e)\n",
            time, spin, mass, actions..., horizon_angular_velocity, imag.(mode_frequencies(solution))...)
    @printf("  %-12s %3s %3s %10s %12s %13s %13s %9s %12s\n", "component q", "m", "l", "ω_c", "ω_c − mΩ_H",
            "F_E into BH", "F_J into BH", "sign ok", "da/dt part")
    total_energy, total_angular_momentum = 0.0, 0.0
    for channel in 1:C.n_channels(mode_set)
        component_index = mode_set.channel_component[channel]
        component = mode_set.components[component_index]
        frequency = component_frequencies[component_index]
        azimuthal = C.channel_m(mode_set, channel)
        wronskian = imag(conj(state.profiles[channel]) * state.slopes[channel])
        outward_energy = prefactor * frequency * wronskian
        outward_angular_momentum = prefactor * azimuthal * wronskian
        total_energy += outward_energy
        total_angular_momentum += outward_angular_momentum
        mass_rate = -coupling_to_hole * outward_energy
        angular_momentum_rate = -coupling_to_hole * mass * outward_angular_momentum
        spin_rate = angular_momentum_rate / mass^2 - 2spin * mass_rate / mass
        detuning = frequency - azimuthal * horizon_angular_velocity
        into_hole = -outward_energy
        sign_ok = into_hole == 0 || sign(into_hole) == sign(frequency * detuning)
        @printf("  %-12s %3d %3d %10.6f %12.4e %13.4e %13.4e %9s %12.4e\n", string(component.q), azimuthal,
                mode_set.channel_l[channel], frequency, detuning, into_hole, -outward_angular_momentum,
                sign_ok ? "yes" : "NO", spin_rate)
    end
    reference = outward_flux(solution, 1)
    @printf("  channel sum − outward_flux(solution, 1): energy %.2e, J %.2e (relative to |F_E| %.2e)\n",
            total_energy - reference.energy, total_angular_momentum - reference.angular_momentum, abs(reference.energy))
    @printf("  hole: M' = %.4e, a' = %.4e per M₀\n", -coupling_to_hole * reference.energy,
            (-coupling_to_hole * mass * reference.angular_momentum) / mass^2 - 2spin * (-coupling_to_hole * reference.energy) / mass)
    flush(stdout)
end
