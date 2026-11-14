# Growth-rate convention for driven components, and the frozen-snapshot action budget.
#
# A component with integer vector q contains conj(φ_i) for every q_i < 0, so its complex frequency is
# Σ_i q_i ω_i + i Σ_i |q_i| γ_i (Omiya's 2ω̃₁ − ω̃₂* for q = (2, −1)), not q·ω̃.
# This script solves interacting snapshots (a = 0.9, α = 0.3, cosine, |211⟩ + |322⟩) and records, per point:
#   - the complex mode frequencies, the component frequencies actually used by the solver,
#   - the wave actions N_i, their losses L_{N_i}, and the per-component charges n_q,
#   - the frozen-snapshot budget: −L_{N_i} against 2γ_i N_i and against Σ_q q_i 2 Im ω̃_q n_q,
#   - the profiles, so that two runs (before and after a code change) can be compared.
#
#   julia --project=. validation/conjugate_growth_check.jl save <output.jls>
#   julia --project=. validation/conjugate_growth_check.jl compare <before.jls> <after.jls>
using AxionCloud, Printf, LinearAlgebra, Serialization
const AC = AxionCloud

const target_amplitudes = [[1e-2, 1e-4], [1e-2, 1e-3]]

"""Charge n_q of every Fourier component between the numerical endpoints (N_i = Σ_q q_i n_q)."""
function component_actions(solution)
    problem, layout = solution.problem, solution.layout
    frequencies = AC.real_component_frequencies(problem, layout, solution.unknowns)
    weights = AC.simpson_weights(problem.mesh.rstar)
    sum(weights[node] .* AC.component_action_densities(problem, problem.mesh.offsets[node],
                                                         AC.physical_state(layout, solution.unknowns, solution.scale, node).profiles,
                                                         frequencies)
        for node in 1:layout.n_nodes)
end

function record(solution)
    mode_set = solution.problem.mode_set
    complex_component_frequencies = AC.component_frequencies(mode_set, mode_frequencies(solution))
    component_n = component_actions(solution)
    growth = growth_rates(solution)
    actions, action_loss = phase_actions(solution), action_losses(solution)
    weights = AC.phase_weights(mode_set)
    frozen_budget = weights * (2 .* imag.(complex_component_frequencies) .* component_n)
    (amplitudes=copy(solution.amplitudes), mode_frequencies=mode_frequencies(solution),
     component_q=[c.q for c in mode_set.components], component_frequencies=complex_component_frequencies,
     component_actions=component_n, actions=actions, action_losses=action_loss,
     two_gamma_N=2 .* growth .* actions, frozen_budget=frozen_budget,
     profiles=profiles(solution), iterations=solution.iterations)
end

function solve_points()
    problem = RadialProblem(KerrBackground(0.9), CosinePotential(), 0.3, ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)]))
    linear = solve_radial(problem, [1e-6, 1e-7])
    records = []
    for target in target_amplitudes
        solution = linear
        for fraction in exp10.(range(-4, 0; length=9)[2:end])
            solution = solve_radial(problem, target .* fraction; guess=solution)
        end
        push!(records, record(solution))
    end
    records
end

function report(record)
    println("== A = ", record.amplitudes, "   Newton iterations = ", record.iterations)
    for (index, mode_frequency) in enumerate(record.mode_frequencies)
        @printf("  mode %d: ω = %.12f  γ = %+.6e\n", index, real(mode_frequency), imag(mode_frequency))
    end
    println("  components (q, Re ω̃_q, Im ω̃_q, n_q):")
    for (q, frequency, n) in zip(record.component_q, record.component_frequencies, record.component_actions)
        @printf("    %-10s %.8f  %+.4e  %.4e\n", string(q), real(frequency), imag(frequency), n)
    end
    for i in eachindex(record.actions)
        @printf("  phase %d: N = %.6e  −L_N = %+.6e  2γN = %+.6e  Σ q 2Im(ω̃_q) n_q = %+.6e\n", i, record.actions[i],
                -record.action_losses[i], record.two_gamma_N[i], record.frozen_budget[i])
        @printf("           (−L_N)/(2γN) − 1 = %+.2e   (−L_N)/(Σ q 2Im ω̃ n) − 1 = %+.2e\n",
                -record.action_losses[i] / record.two_gamma_N[i] - 1, -record.action_losses[i] / record.frozen_budget[i] - 1)
    end
end

relative(a, b) = norm(a - b) / norm(b)

function compare(before, after)
    for (old, new) in zip(before, after)
        println("== A = ", new.amplitudes)
        for i in eachindex(new.mode_frequencies)
            @printf("  mode %d: Δω/ω = %+.2e   Δγ/γ = %+.2e\n", i,
                    real(new.mode_frequencies[i] - old.mode_frequencies[i]) / real(old.mode_frequencies[i]),
                    imag(new.mode_frequencies[i] - old.mode_frequencies[i]) / imag(old.mode_frequencies[i]))
        end
        @printf("  profiles, all channels: ‖Δ‖/‖·‖ = %.2e\n", relative(new.profiles, old.profiles))
        for channel in axes(new.profiles, 2)
            @printf("    channel %2d: ‖Δ‖/‖·‖ = %.2e   (‖R‖ = %.2e)\n", channel,
                    relative(new.profiles[:, channel], old.profiles[:, channel]), norm(old.profiles[:, channel]))
        end
        for i in eachindex(new.actions)
            @printf("  phase %d: ΔN/N = %+.2e   ΔL_N/L_N = %+.2e\n", i, (new.actions[i] - old.actions[i]) / old.actions[i],
                    (new.action_losses[i] - old.action_losses[i]) / old.action_losses[i])
        end
    end
end

if ARGS[1] == "save"
    records = solve_points()
    foreach(report, records)
    serialize(ARGS[2], records)
elseif ARGS[1] == "compare"
    before, after = deserialize(ARGS[2]), deserialize(ARGS[3])
    println("before:"); foreach(report, before)
    println("\nafter:"); foreach(report, after)
    println("\nafter − before:"); compare(before, after)
end
