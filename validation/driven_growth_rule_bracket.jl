# Bracket the prescription for the imaginary part of driven components in terms of the amplitude rates.
#
# A driven Fourier component with integer vector q has no eigenvalue of its own; its complex frequency
# is a prescription. Three rules are available (`driven_growth_rule` of `RadialProblem`):
#   :conjugate  Im ω̃_q = Σ_i |q_i| γ_i   (default; Omiya's 2ω̃₁ − ω̃₂* for q = (2, −1))
#   :signed     Im ω̃_q = Σ_i  q_i  γ_i   (ω̃_q = q·ω̃, the rule before commit 41fc2e7)
#   :zero       Im ω̃_q = 0                (real frequency q·ω)
# Fundamental components keep their own eigenvalues under all three.
#
# For each point this solves the snapshot under each rule and prints ω_i, γ_i, the wave actions N_i
# (`phase_actions`), their losses L_{N_i} (`action_losses`) and the fixed-background action-balance
# rates Ȧ_i (`amplitude_rates(solution; balance=:action)`), then the relative differences of each rule
# from :conjugate. Everything is in units of the current black-hole mass M (ω = Mω̂, rates per t/M).
#
# Points:
#   P1-P4  a = 0.9, α = 0.3, cosine, A = (1e-2, 1e-4) and (1e-2, 1e-3): the two points of
#          validation/conjugate_growth_check.jl, on two meshes: that script's default
#          `RadialMesh(background)` and the contained mesh `RadialMesh(background, α, mode_set)`
#          of the spindown driver.
#   P5-P7  points read from the stored cosine spindown runs, on the spindown driver's mesh. The stored
#          amplitudes A are the solver's own amplitudes (scripts/evolution_spindown.jl passes them
#          straight to `solve_radial`; f_a enters only the black-hole coupling), so no conversion is
#          needed. The radial problem is rebuilt as `evolve_with_spindown` builds it at time τ:
#          `with_background(problem_initial, KerrBackground(a), α₀ b)` with `problem_initial` made for the
#          run's initial spin on `RadialMesh(KerrBackground(a₀), α₀, mode_set)`, and (a, b) read from the
#          stored row. The fixed-background action balance is evaluated there, without the spindown
#          drift terms.
#
#   julia --project=. validation/driven_growth_rule_bracket.jl [point numbers, default all]
using AxionCloud, Printf, DelimitedFiles
const AC = AxionCloud

const rules = (:conjugate, :signed, :zero)
const mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
const results_directory = joinpath(@__DIR__, "..", "results")

"""Row of a stored spindown diagnostics CSV chosen by `selector` (:max_A1, :max_A2, :max_ratio = A2/A1)."""
function stored_row(file_name, selector)
    table, header = readdlm(joinpath(results_directory, file_name), ','; header=true)
    column(name) = table[:, findfirst(==(name), vec(header))]
    keep = findall(i -> isfinite(column("A1")[i]) && isfinite(column("A2")[i]), axes(table, 1))
    scores = selector == :max_A1 ? column("A1")[keep] : selector == :max_A2 ? column("A2")[keep] :
             column("A2")[keep] ./ column("A1")[keep]
    row = keep[argmax(scores)]
    (time=column("time")[row], amplitudes=[column("A1")[row], column("A2")[row]], mass=column("b")[row], spin=column("a")[row])
end

"""Amplitude continuation from 1e-3 A to A in 10 logarithmic steps (the spindown restart uses 30), for one rule; returns the last solution."""
function solve_with_continuation(problem, amplitudes)
    solution = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=10))
        solution = solve_radial(problem, amplitudes .* factor; guess=solution)
        solution.converged || return solution
    end
    solution
end

"""
Solve under one rule. With `warm_start` (the converged solution of another rule at the same amplitudes, same
unknown layout) Newton starts from it directly; if that does not converge, amplitude continuation starts from it
at 1e-3 A. Without it, continuation starts cold. The converged solution is the same fixed point either way.
"""
function solve_under_rule(problem, amplitudes, warm_start)
    warm_start === nothing && return solve_with_continuation(problem, amplitudes), "continuation from cold"
    solution = solve_radial(problem, amplitudes; guess=warm_start)
    solution.converged && return solution, "warm start from conjugate solution"
    return solve_with_continuation(problem, amplitudes), "warm start failed; continuation from cold"
end

struct Point
    label::String
    amplitudes::Vector{Float64}
    problem_builder::Function   # rule -> RadialProblem
    note::String
end

function fixed_spin_point(label, amplitudes, mesh_kind)
    builder = function (rule)
        background = KerrBackground(0.9)
        mesh = mesh_kind == :default ? RadialMesh(background) : RadialMesh(background, 0.3, mode_set)
        RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, driven_growth_rule=rule)
    end
    mesh_text = mesh_kind == :default ? "default mesh RadialMesh(background)" : "contained mesh RadialMesh(background, α, mode_set)"
    Point(label, amplitudes, builder, "a = 0.9, α = 0.3, cosine, $mesh_text")
end

function stored_point(label, file_name, selector, initial_spin)
    row = stored_row(file_name, selector)
    builder = function (rule)
        background0 = KerrBackground(initial_spin)
        mesh = RadialMesh(background0, 0.3, mode_set)
        initial = RadialProblem(background0, CosinePotential(), 0.3, mode_set; mesh, driven_growth_rule=rule)
        with_background(initial, KerrBackground(row.spin), 0.3 * row.mass)
    end
    Point(label, row.amplitudes, builder,
          @sprintf("from %s (%s): τ = %.6e, a = %.6f, b = %.6f, α = %.6f, cosine; spindown-driver mesh of a₀ = %g",
                   file_name, selector, row.time, row.spin, row.mass, 0.3 * row.mass, initial_spin))
end

function point_list()
    [fixed_spin_point("P1 A=(1e-2,1e-4), default mesh", [1e-2, 1e-4], :default),
     fixed_spin_point("P2 A=(1e-2,1e-3), default mesh", [1e-2, 1e-3], :default),
     fixed_spin_point("P3 A=(1e-2,1e-4), contained mesh", [1e-2, 1e-4], :contained),
     fixed_spin_point("P4 A=(1e-2,1e-3), contained mesh", [1e-2, 1e-3], :contained),
     stored_point("P5 largest A1, spindown a0=0.99 fa=3e-3", "spindown_a0.99_alpha0.3_fa3e-3_diagnostics.csv", :max_A1, 0.99),
     stored_point("P6 largest A2/A1, spindown a0=0.99 fa=1e-2", "spindown_a0.99_alpha0.3_fa1e-2_diagnostics.csv", :max_ratio, 0.99),
     stored_point("P7 largest A2, spindown a0=0.99 fa=3e-3", "spindown_a0.99_alpha0.3_fa3e-3_diagnostics.csv", :max_A2, 0.99)]
end

relative_difference(value, reference) = (value - reference) / reference

function run_point(point)
    println("\n", "="^100, "\n", point.label, "\n  ", point.note); flush(stdout)
    @printf("  A = (%.10e, %.10e)\n", point.amplitudes...)
    records = Dict{Symbol, Any}()
    conjugate_solution = nothing
    for rule in rules
        problem = point.problem_builder(rule)
        solve_method = ""
        warm_start = haskey(records, :conjugate) ? conjugate_solution : nothing
        elapsed = @elapsed solution, solve_method = solve_under_rule(problem, point.amplitudes, warm_start)
        rule == :conjugate && (conjugate_solution = solution)
        if !solution.converged
            @printf("  rule %-9s: NOT CONVERGED (|F| = %.3e after %d iterations); skipped\n", rule, solution.residual_norm, solution.iterations)
            continue
        end
        frequencies = mode_frequencies(solution)
        actions, losses = phase_actions(solution), action_losses(solution)
        rates = amplitude_rates(solution; balance=:action)
        component_index = findfirst(component -> component.q == [-1, 2], problem.mode_set.components)
        driven_frequency = AC.component_frequencies(problem.mode_set, frequencies, rule)[component_index]
        threshold_gap = real(driven_frequency) - problem.mass_coupling
        records[rule] = (frequencies=frequencies, actions=actions, losses=losses, rates=rates)
        @printf("  rule %-9s (%.0f s; %s; Newton iterations in last solve %d)\n", rule, elapsed, solve_method, solution.iterations)
        @printf("    q=(-1,2): Re ω̃_q − α = %.6e, Im ω̃_q = %+.4e, ratio Im/gap = %+.3e\n", threshold_gap, imag(driven_frequency),
                imag(driven_frequency) / threshold_gap)
        for i in 1:2
            @printf("    mode %d: ω = %.12f  γ = %+.6e  N = %.6e  L_N = %+.6e  Ȧ = %+.6e  (Ȧ/A = %+.6e)\n", i, real(frequencies[i]),
                    imag(frequencies[i]), actions[i], losses[i], rates[i], rates[i] / point.amplitudes[i])
        end
        flush(stdout)
    end
    haskey(records, :conjugate) || return nothing
    reference = records[:conjugate]
    for rule in (:signed, :zero)
        haskey(records, rule) || continue
        record = records[rule]
        println("  relative difference ", rule, " vs conjugate, (X_rule − X_conj)/X_conj:")
        for i in 1:2
            @printf("    mode %d: Δω/ω = %+.2e   Δγ/γ = %+.2e   ΔN/N = %+.2e   ΔL_N/L_N = %+.2e   ΔȦ/Ȧ = %+.2e\n", i,
                    relative_difference(real(record.frequencies[i]), real(reference.frequencies[i])),
                    relative_difference(imag(record.frequencies[i]), imag(reference.frequencies[i])),
                    relative_difference(record.actions[i], reference.actions[i]),
                    relative_difference(record.losses[i], reference.losses[i]),
                    relative_difference(record.rates[i], reference.rates[i]))
        end
    end
    return records
end

selected = isempty(ARGS) ? nothing : parse.(Int, ARGS)
for (index, point) in enumerate(point_list())
    (selected === nothing || index in selected) && run_point(point)
end
