# Size of the nonlinear source N_{l,q} at the numerical endpoints, where the boundary data are the *linear*
# asymptotic expansions (outline, "Linear boundary conditions"). For every channel the local measure is
#     |N| / |B R|,
# the neglected source against the linear term it competes with in (Δ R')' + B R + a²(ω²−α²) C R + N = 0
# (both multiplied by Δ, as in `radial_rhs`, which is finite at the horizon). The ratio is also printed a few
# nodes inside each endpoint to show how fast it falls off.
# Run from the package root:  julia --project=. validation/nonlinear_boundary_terms.jl
using AxionCloud, Printf
const AC = AxionCloud

"""Per-channel |ΔN| / |ΔB R| at mesh `node` of a converged solution."""
function nonlinear_to_linear_ratios(solution, node)
    problem, layout = solution.problem, solution.layout
    background, mode_set, mass_coupling = problem.background, problem.mode_set, problem.mass_coupling
    a = background.a
    s = problem.mesh.offsets[node]
    r = AC.radius(background, s)
    horizon_function, sum_r2a2 = AC.delta(background, s), AC.r2a2(background, s)
    state = AC.physical_state(layout, solution.unknowns, solution.scale, node)
    source = AC.nonlinear_source(problem.quadrature, problem.potential, background, mass_coupling, s, state.profiles)
    frequencies = AC.component_frequencies(mode_set, mode_frequencies(solution))
    ratios = zeros(AC.n_channels(mode_set))
    for (component_index, component) in enumerate(mode_set.components)
        frequency = frequencies[component_index]
        for (local_index, channel) in enumerate(mode_set.component_channels[component_index])
            l = component.ls[local_index]
            delta_B = (sum_r2a2 * frequency - a * component.m)^2 +
                      horizon_function * (2a * component.m * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - l * (l + 1))
            ratios[channel] = abs(horizon_function * source[channel]) / abs(delta_B * state.profiles[channel])
        end
    end
    return ratios
end

function report(label, solution)
    problem, layout = solution.problem, solution.layout
    mode_set, mesh = problem.mode_set, problem.mesh
    println("\n== ", label, ": A = ", solution.amplitudes, ", ω = ", mode_frequencies(solution),
            ", r* ∈ [", mesh.rstar[1], ", ", mesh.rstar[end], "]")
    n = layout.n_nodes
    sample_nodes = (("horizon", 1), ("horizon + 25 nodes", 26), ("outer − 50 nodes", n - 50), ("outer − 25 nodes", n - 25), ("outer", n))
    columns = [nonlinear_to_linear_ratios(solution, node) for (_, node) in sample_nodes]
    @printf("  %-22s %-5s", "component q", "l")
    for (name, node) in sample_nodes
        @printf(" %12s", name)
    end
    println()
    @printf("  %-22s %-5s", "", "")
    for (_, node) in sample_nodes
        @printf("   r*=%7.1f", mesh.rstar[node])
    end
    println()
    for channel in 1:AC.n_channels(mode_set)
        component = AC.channel_component(mode_set, channel)
        @printf("  %-22s %-5d", string(component.q), mode_set.channel_l[channel])
        for column in columns
            @printf(" %12.1e", column[channel])
        end
        println()
    end
    @printf("  largest at the horizon endpoint %.1e, at the outer endpoint %.1e\n", maximum(columns[1]), maximum(columns[end]))
    flush(stdout)
end

"""Geometric continuation in the amplitudes from a linear-regime solution to `target`."""
function continue_to(problem, target; steps=10)
    solution = solve_radial(problem, target .* 1e-4)
    for fraction in exp10.(range(-4, 0; length=steps + 1)[2:end])
        solution = solve_radial(problem, target .* fraction; guess=solution)
    end
    return solution
end

mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])

# 1. a = 0.99, α = 0.3, E₁ = E₂ = 650 (joint solution of the profile comparison), r* ∈ [−100, 200], step 0.4.
background = KerrBackground(0.99)
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=200.0, step=0.4)
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
joint = amplitudes_for_energy(problem, [650.0, 650.0])
report("a = 0.99, E₁ = E₂ = 650 joint solution", joint)

# 2. a = 0.9, α = 0.3, default mesh: states of the long fixed-background evolution.
problem = RadialProblem(KerrBackground(0.9), CosinePotential(), 0.3, mode_set)
report("a = 0.9, peak of |211⟩ in the long evolution (t = 6.58e8)", continue_to(problem, [1.905e-2, 9.6e-6]))
report("a = 0.9, saturated state of the long evolution (t = 1.5e9)", continue_to(problem, [5.52e-3, 3.44e-5]))
