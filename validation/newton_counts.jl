# How many residuals, Jacobian assemblies and factorisations does each radial solve cost along an evolution?
# Counts, not seconds, are the metric (the cost of one residual or Jacobian changes with the assembly code).
#
# The counts come from wrappers, defined in this script only, around `residual`, `jacobian`, `ScaledFactorization`
# and `\` (more specific methods that record an event and `invoke` the library method). Per solve the event string
# reads, e.g., `RSRSRJF`: R = residual, J = Jacobian assembly, F = factorisation, S = back-substitution for a step.
# Derived columns: newton = steps taken with a fresh factorisation ("FS"), chord = the other steps, back = extra
# residuals of the line search, finJ = whether the solve ended with the refactorisation at the converged point.
#
# Regimes: four states of the f_a/M_Pl = 0.0466873 spindown run (log
# results/spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart-61110567.out): linear (seed), |211⟩ plateau (t = 1.5e9),
# transition (t = 3.1e10, where the log shows 12-13 iterations per solve) and both clouds nonlinear (t = 1e11).
# Each state is built by amplitude continuation on its background. For each state the script runs a short
# evolve_with_spindown segment (`steps` accepted steps of Tsit5) with the tangent predictor and without it.
#
# Run from the package root (about 30 min; submit to Slurm):
#   julia --project=. validation/newton_counts.jl [states, default linear,plateau,transition,nonlinear] [steps, default 10]
using AxionCloud, LinearAlgebra, SparseArrays, Printf
const AC = AxionCloud

# --- event recording ------------------------------------------------------------------------------------------
const RECORDING = Ref(false)
const EVENTS = Char[]
record!(event::Char) = (RECORDING[] && push!(EVENTS, event); nothing)

function AC.residual(problem::AC.RadialProblem, layout::AC.UnknownLayout, unknowns::Vector{Float64},
                     amplitudes::Vector{Float64}, scale::Float64)
    record!('R')
    invoke(AC.residual, Tuple{AC.RadialProblem, AC.UnknownLayout, AbstractVector, AbstractVector, Any},
           problem, layout, unknowns, amplitudes, scale)
end
function AC.jacobian(problem::AC.RadialProblem, layout::AC.UnknownLayout, unknowns::Vector{Float64},
                     amplitudes::Vector{Float64}, scale::Float64)
    record!('J')
    invoke(AC.jacobian, Tuple{AC.RadialProblem, AC.UnknownLayout, AbstractVector, AbstractVector, Any},
           problem, layout, unknowns, amplitudes, scale)
end
function AC.ScaledFactorization(matrix::SparseMatrixCSC{Float64, Int}, scaling::AC.NewtonScaling)
    record!('F')
    invoke(AC.ScaledFactorization, Tuple{AbstractMatrix, AC.NewtonScaling}, matrix, scaling)
end
function Base.:\(scaled::AC.ScaledFactorization, right_hand_side::Vector{Float64})
    record!('S')
    invoke(Base.:\, Tuple{AC.ScaledFactorization, AbstractVecOrMat}, scaled, right_hand_side)
end

"""Events, iterations and outcome of one recorded solve."""
struct SolveRecord
    events::String
    iterations::Int
    converged::Bool
end
const SOLVE_LOG = SolveRecord[]

function AC.solve_radial(problem::AC.RadialProblem, amplitudes::Vector{Float64}; kwargs...)
    # nested solves (single-mode frequency solves of a cold start) are counted in the outer record
    RECORDING[] && return invoke(AC.solve_radial, Tuple{AC.RadialProblem, AbstractVector}, problem, amplitudes; kwargs...)
    empty!(EVENTS)
    RECORDING[] = true
    solution = try
        invoke(AC.solve_radial, Tuple{AC.RadialProblem, AbstractVector}, problem, amplitudes; kwargs...)
    finally
        RECORDING[] = false
    end
    push!(SOLVE_LOG, SolveRecord(String(copy(EVENTS)), solution.iterations, solution.converged))
    return solution
end

"""Counts derived from an event string (see the header)."""
function event_counts(events::AbstractString)
    newton_steps = length(collect(eachmatch(r"FS", events)))
    trials = [length(m.match) for m in eachmatch(r"(?<=S)R+", events)]
    (residuals=count(==('R'), events), jacobians=count(==('J'), events), factorizations=count(==('F'), events),
     newton_steps, chord_steps=count(==('S'), events) - newton_steps, backtracks=sum(trials .- 1; init=0),
     final_refactorization=endswith(events, "JF"))
end

function print_summary(label, records)
    counts = [event_counts(record.events) for record in records]
    per_solve(field) = sum(getfield.(counts, field); init=0) / length(records)
    @printf("%-34s solves %3d | per solve: iterations %5.2f  R %5.2f  J %4.2f  F %4.2f  newton %4.2f  chord %5.2f  back %5.2f  finJ %4.2f | failed %d\n",
            label, length(records), sum(record.iterations for record in records) / length(records), per_solve(:residuals),
            per_solve(:jacobians), per_solve(:factorizations), per_solve(:newton_steps), per_solve(:chord_steps),
            per_solve(:backtracks), count(count -> count.final_refactorization, counts) / length(records),
            count(record -> !record.converged, records))
    pattern_counts = Dict{String, Int}()
    for record in records
        pattern_counts[record.events] = get(pattern_counts, record.events, 0) + 1
    end
    for (pattern, occurrences) in first(sort(collect(pattern_counts); by=last, rev=true), 5)
        @printf("    %3d × %s\n", occurrences, length(pattern) > 100 ? pattern[1:100] * "…" : pattern)
    end
    flush(stdout)
end

# --- states -------------------------------------------------------------------------------------------------------
const DECAY_CONSTANT_RATIO = 0.0466873
const MODE_SET = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
const REFERENCE = let background = KerrBackground(0.9)
    RadialProblem(background, CosinePotential(), 0.3, MODE_SET; mesh=RadialMesh(background, 0.3, MODE_SET))
end
# (amplitudes, b, j, t) from the log
const STATES = Dict(
    "linear" => ([9.85578501492328e-7, 3.057454189623335e-9], 1.0, 0.9, 0.0),
    "plateau" => ([0.00280773553, 3.3961332e-9], 0.9814753887857007, 0.8374944341429903, 1.5e9),
    "transition" => ([0.0024024134667967956, 2.4071466062915234e-6], 0.9854426851242335, 0.8458838283206854,
                     3.0968105831506714e10),
    "nonlinear" => ([1.4719255751539253e-4, 4.100281105490423e-6], 0.9954602054965745, 0.8671482601144008, 1e11),
)
problem_at(state::BackgroundState) =
    with_background(REFERENCE, KerrBackground(state.angular_momentum / state.mass^2), REFERENCE.mass_coupling * state.mass)

"""Converged solution at a state by amplitude continuation on its background (cold for the linear state)."""
function state_solution(amplitudes, background)
    problem = problem_at(background)
    maximum(amplitudes) < 1e-5 && return solve_radial(problem, amplitudes)
    solution = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=30))
        solution = solve_radial(problem, amplitudes .* factor; guess=solution)
        solution.converged || error("continuation failed at factor $factor")
    end
    return solution
end

state_names = length(ARGS) ≥ 1 ? split(ARGS[1], ',') : ["linear", "plateau", "transition", "nonlinear"]
steps = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 10
for name in state_names
    amplitudes, mass, angular_momentum, time = STATES[name]
    background = BackgroundState(mass, angular_momentum)
    solution = state_solution(amplitudes, background)
    println("\n=== $name: t = $time, A = $amplitudes, b = $mass, j = $angular_momentum")
    for tangent_predictor in (true, false)
        start = length(SOLVE_LOG) + 1
        result = evolve_with_spindown(REFERENCE, amplitudes, (time, time + 1e13), DECAY_CONSTANT_RATIO;
                                      initial_background=background, initial_guess=solution, maxiters=steps, tangent_predictor)
        stats = result.ode_solution.stats
        @printf("tangent_predictor = %s: %d accepted, %d rejected steps, t = %.6e\n", tangent_predictor, stats.naccept,
                stats.nreject, result.ode_solution.t[end])
        # the first two solves are at the starting state (initial solve and the integrator's first evaluation)
        print_summary("$name, tangent_predictor = $tangent_predictor", SOLVE_LOG[start+2:end])
    end
end
