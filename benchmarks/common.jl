"""
Shared setup of the benchmark and regression harness: the reference states, their
stored states, and the thin wrappers around the library internals that
the harness calls. `regression.jl`, `timing.jl` and `noise_floor.jl` include this
file and nothing else from the package's own scripts.

The library functions used that are not exported from `AxionCloud` are listed in
`INTERNAL_INTERFACE` below; if an optimisation changes one of their signatures,
the wrapper in this file is the only place to adapt.

Notation introduced here (not used elsewhere in the repo):

* "reference problem": the `RadialProblem` at the initial background (`b = 1`,
  `j = a₀`) that defines the units and the mesh of an evolution. The evolution
  drivers call it `problem`.
* "state": a converged radial solution at given amplitudes and black-hole
  `BackgroundState(b, j)`, together with the reference problem it belongs to.
"""
module BenchmarkHarness

using AxionCloud
using LinearAlgebra, Printf, Dates, SparseArrays

const AC = AxionCloud

const BENCHMARK_DIRECTORY = @__DIR__
const DECAY_CONSTANT_RATIO = 0.0466873

# Library internals reached through `AC.` (not exported):
#   residual, jacobian, UnknownLayout, NewtonScaling, ScaledFactorization, RadialParameters,
#   parameter_range, frequency_range, n_unknowns, _complexify, _spindown_rhs!, current_problem,
#   EvolutionCache fields (accepted_solution, last_solution, diagnostics)
const INTERNAL_INTERFACE = (:residual, :jacobian, :UnknownLayout, :NewtonScaling, :ScaledFactorization,
                            :RadialParameters, :parameter_range, :frequency_range, :n_unknowns, :_complexify,
                            :_spindown_rhs!, :current_problem)

# --- state definitions ---------------------------------------------------------

"""
    StateSpec

Definition of one reference state. `construction` is `:cold` (one cold solve at
`amplitudes`) or `:continuation` (cold solve at `1e-3 × amplitudes` on the target
background, then geometric amplitude continuation to `amplitudes`, as the restart
branch of that script does). `solve_tolerance` is the Newton tolerance used to
build the state and for the perturbed solves on it.
"""
struct StateSpec
    name::String
    description::String
    mass_coupling::Float64
    initial_spin::Float64
    amplitudes::Vector{Float64}
    black_hole_mass::Float64
    black_hole_angular_momentum::Float64
    construction::Symbol
    solve_tolerance::Float64
    run_tolerance::Float64   # Newton tolerance of the perturbed and same-point solves (1e-10 is the evolution default; unreachable on the 1355-node mesh)
    include_evolution::Bool
end

background_state(spec::StateSpec) = BackgroundState(spec.black_hole_mass, spec.black_hole_angular_momentum)

const STATE_SPECS = [
    StateSpec("linear", "linear regime, the seed amplitudes printed by scripts/evolution_spindown.jl (alpha 0.3, a 0.9)",
              0.3, 0.9, [9.85578501492328e-7, 3.057454189623335e-9], 1.0, 0.9, :cold, 1e-12, 1e-10, false),
    StateSpec("plateau", "|211> plateau of the f_a/M_Pl = 0.0466873 run at t = 1.5e9 (log 61110567, first line)",
              0.3, 0.9, [0.00280773553, 3.3961332e-9], 0.9814753887857007, 0.8374944341429903, :continuation, 1e-12, 1e-10, true),
    StateSpec("nonlinear", "both clouds nonlinear, t = 1e11 of the same log (last line)",
              0.3, 0.9, [1.4719255751539253e-4, 4.100281105490423e-6], 0.9954602054965745, 0.8671482601144008,
              :continuation, 1e-12, 1e-10, false),
    StateSpec("large_mesh", "linear regime at alpha 0.187, a 0.9 (1355 mesh nodes), the mesh size production is heading to",
              0.187, 0.9, [1e-6, 1e-6], 1.0, 0.9, :cold, 1e-9, 1e-8, false),
]

spec_named(name::AbstractString) = STATE_SPECS[findfirst(spec -> spec.name == name, STATE_SPECS)]

# --- problem setup -------------------------------------------------------------

const MODE_SET = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])

"""Reference problem (initial background, production mesh, cosine potential) of a state."""
function reference_problem(spec::StateSpec)
    background = KerrBackground(spec.initial_spin)
    mesh = RadialMesh(background, spec.mass_coupling, MODE_SET)
    return RadialProblem(background, CosinePotential(), spec.mass_coupling, MODE_SET; mesh)
end

"""The reference problem on the background `(b, j)`, in units of the current mass (what `evolve_with_spindown` solves)."""
problem_at(reference::RadialProblem, state::BackgroundState) =
    with_background(reference, KerrBackground(state.angular_momentum / state.mass^2), reference.mass_coupling * state.mass)

"""
    BenchmarkState

`spec`, the reference problem, the problem at the state's background, the
background, and the converged `RadialSolution` there.
"""
struct BenchmarkState
    spec::StateSpec
    reference::RadialProblem
    problem::RadialProblem
    background::BackgroundState{Float64}
    solution::RadialSolution
end

"""
    solution_from_unknowns(problem, amplitudes, unknowns, scale)

The `RadialSolution` at a given unknown vector, with a fresh factorisation of the
Jacobian at exactly that point (`solve_radial` may return one taken within `√tol`
of it). This is the one place where a stored state becomes a library object; adapt it if
the solution type or its factorisation changes.
"""
function solution_from_unknowns(problem::RadialProblem, amplitudes::AbstractVector, unknowns::AbstractVector, scale)
    layout = AC.UnknownLayout(problem)
    frequencies = AC._complexify(unknowns[AC.parameter_range(layout)][AC.frequency_range(layout)])
    scaling = AC.NewtonScaling(problem, layout, frequencies)
    factorization = AC.ScaledFactorization(AC.jacobian(problem, layout, unknowns, amplitudes, scale), scaling)
    scaled_residual = norm(scaling.row .* AC.residual(problem, layout, unknowns, amplitudes, scale), Inf)
    parameters = AC.RadialParameters(layout, unknowns[AC.parameter_range(layout)], amplitudes, scale)
    return RadialSolution(problem, layout, collect(float.(amplitudes)), float(scale), collect(float.(unknowns)), parameters,
                          factorization, 0, scaled_residual, true)
end

# --- state construction --------------------------------------------------------

"""Solve to `tolerance`, failing loudly if the solve does not converge."""
function converged_solve(problem::RadialProblem, amplitudes; guess=nothing, tolerance)
    solution = solve_radial(problem, amplitudes; guess, tol=tolerance)
    solution.converged || error("solve_radial did not converge at amplitudes $amplitudes " *
                                "(|F| = $(solution.residual_norm), tolerance $tolerance)")
    return solution
end

"""Cold solve at `1e-3 × amplitudes`, then `steps` geometrically spaced amplitude factors up to `1`."""
function continuation_solution(problem::RadialProblem, amplitudes, tolerance; steps=30, verbose=false)
    guess = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=steps))
        guess = converged_solve(problem, amplitudes .* factor; guess, tolerance)
        verbose && @printf("    continuation factor %.4e: %d iterations, |F| = %.2e\n", factor, guess.iterations, guess.residual_norm)
    end
    return guess
end

"""
    build_state(spec; continuation_steps=30, stage_tolerance=max(spec.solve_tolerance, 1e-10), verbose=false)

Build a state from scratch. The construction runs at `stage_tolerance`:
a cold solve cannot go below a scaled residual of about `3e-11`, because the Newton
scaling is built from the starting frequencies. A final warm solve at
`spec.solve_tolerance` (scaling now built from the converged frequencies, residual
near its roundoff floor) makes rebuilding on another machine change the unknowns
only at roundoff level; `noise_floor.jl` measures that.
"""
function build_state(spec::StateSpec; continuation_steps=30, stage_tolerance=max(spec.solve_tolerance, 1e-10), verbose=false)
    reference = reference_problem(spec)
    state = background_state(spec)
    problem = problem_at(reference, state)
    raw = if spec.construction == :cold
        converged_solve(problem, spec.amplitudes; tolerance=stage_tolerance)
    elseif spec.construction == :continuation
        continuation_solution(problem, spec.amplitudes, stage_tolerance; steps=continuation_steps, verbose)
    else
        error("unknown construction $(spec.construction)")
    end
    polished = converged_solve(problem, raw.amplitudes; guess=raw, tolerance=spec.solve_tolerance)
    solution = solution_from_unknowns(problem, polished.amplitudes, polished.unknowns, polished.scale)
    return BenchmarkState(spec, reference, problem, state, solution)
end

# --- stored states ---------------------------------------------------------------

"""
Each state is stored as the raw little-endian `Float64` unknown vector in
`benchmarks/states/<name>.bin`, tracked in git. Everything else that defines the
state (amplitudes, `b`, `j`, mass coupling, mesh) is in the `StateSpec`, and the
unknown scale is `maximum(amplitudes)`. These files are fixed inputs: the
regression references were computed from them, and rebuilding them from scratch
(`build_states.jl`) gives a different state at the level of the Newton tolerance,
which is not comparable with the references.
"""
states_directory() = get(ENV, "BENCH_STATES_DIR", joinpath(BENCHMARK_DIRECTORY, "states"))
state_path(spec::StateSpec) = joinpath(states_directory(), "$(spec.name).bin")

function save_state(state::BenchmarkState)
    mkpath(states_directory())
    write(state_path(state.spec), state.solution.unknowns)
end

"""The stored state, with a `RadialSolution` rebuilt from its unknown vector by `solution_from_unknowns`."""
function load_state(spec::StateSpec)
    path = state_path(spec)
    isfile(path) || error("no stored state at $path; run benchmarks/build_states.jl")
    unknowns = reinterpret(Float64, read(path))
    reference = reference_problem(spec)
    state = background_state(spec)
    problem = problem_at(reference, state)
    length(unknowns) == AC.n_unknowns(AC.UnknownLayout(problem)) ||
        error("$path holds $(length(unknowns)) unknowns but the problem has $(AC.n_unknowns(AC.UnknownLayout(problem)))")
    solution = solution_from_unknowns(problem, spec.amplitudes, collect(unknowns), maximum(spec.amplitudes))
    return BenchmarkState(spec, reference, problem, state, solution)
end

selected_specs(names) = isempty(names) ? STATE_SPECS : [spec_named(name) for name in names]

# --- fixed perturbations and vectors -------------------------------------------

"""Amplitude factors of the fixed small perturbation. The first amplitude is left alone: changing the largest amplitude changes the unknown scale and raises the Newton residual floor on the linear state to about 6e-11, close to the default tolerance 1e-10."""
const PERTURBATION_AMPLITUDE_FACTORS = [1.0, 1.01]
"""Relative changes of `b` and `j` of the fixed perturbation (`a = j/b²` changes by `+5e-7`). Shifts of 1e-5 make the warm solve diverge on the linear state; evolution steps are smaller."""
const PERTURBATION_BACKGROUND_FACTORS = (1 - 1e-6, 1 - 1.5e-6)

perturbed_amplitudes(state::BenchmarkState) = state.solution.amplitudes .* PERTURBATION_AMPLITUDE_FACTORS
perturbed_background(state::BenchmarkState) =
    BackgroundState(state.background.mass * PERTURBATION_BACKGROUND_FACTORS[1],
                    state.background.angular_momentum * PERTURBATION_BACKGROUND_FACTORS[2])

"""The solve used by the regression and timing scripts: perturbed amplitudes and background, warm-started from the state."""
function perturbed_solve(state::BenchmarkState; tolerance=state.spec.run_tolerance, guess=state.solution)
    problem = problem_at(state.reference, perturbed_background(state))
    return solve_radial(problem, perturbed_amplitudes(state); guess, tol=tolerance)
end

"""Same-point solve: the state's own amplitudes, warm-started from the state."""
same_point_solve(state::BenchmarkState; tolerance=state.spec.run_tolerance) =
    solve_radial(state.problem, state.solution.amplitudes; guess=state.solution, tol=tolerance)

"""Deterministic unit-size vector of length `count` (no random numbers)."""
fixed_vector(count::Int; frequency=0.37) = [sin(frequency * index) + 0.25cos(0.011frequency * index^1.1) for index in 1:count]

"""Unknown vector displaced from the converged one by a fixed smooth relative pattern (`1e-3`), for evaluating residual and Jacobian away from the root."""
displaced_unknowns(solution::RadialSolution) =
    solution.unknowns .* (1 .+ 1e-3 .* fixed_vector(length(solution.unknowns); frequency=0.53))

"""Indices spread over the unknown vector (first and last included), used for sampled entries."""
sample_indices(count::Int, samples::Int) = unique(round.(Int, range(1, count; length=samples)))

# --- operations shared by regression and timing --------------------------------

residual_at(state::BenchmarkState, unknowns) =
    AC.residual(state.problem, state.solution.layout, unknowns, state.solution.amplitudes, state.solution.scale)

jacobian_at(state::BenchmarkState, unknowns) =
    AC.jacobian(state.problem, state.solution.layout, unknowns, state.solution.amplitudes, state.solution.scale)

newton_scaling(state::BenchmarkState) =
    AC.NewtonScaling(state.problem, state.solution.layout, mode_frequencies(state.solution))

factorize(matrix, scaling) = AC.ScaledFactorization(matrix, scaling)

"""Fresh evolution cache on the reference problem, as `scripts/evolution_spindown.jl` would build."""
evolution_cache(state::BenchmarkState; newton_tolerance=state.spec.run_tolerance) =
    EvolutionCache(state.reference; decay_constant_ratio=DECAY_CONSTANT_RATIO, newton_tol=newton_tolerance,
                   record_metric_work=true)   # the compared balance terms include the energy pieces

"""One `_spindown_rhs!` call at the state, warm-started from the state's own solution."""
function spindown_rhs_call(state::BenchmarkState, cache=evolution_cache(state))
    mode_count = length(state.solution.amplitudes)
    state_vector = vcat(log.(state.solution.amplitudes), state.background.mass, state.background.angular_momentum)
    cache.accepted_solution = state.solution
    cache.failed_attempt = -1   # a failed call would otherwise make every later call return at once
    rates = zeros(mode_count + 2)
    AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    cache.failed_attempt == -1 || error("spindown rhs call: the radial solve did not converge on state '$(state.spec.name)'")
    empty!(cache.diagnostics)
    return rates
end

"""
    evolution_segment(state; time_span=(1.5e9, 1.9e9))

`evolve_with_spindown` on the reference problem from the state over `time_span`;
returns the result and the right-hand-side call count of the integrator.
"""
function evolution_segment(state::BenchmarkState; time_span=(1.5e9, 1.9e9), newton_tolerance=1e-10)
    result = evolve_with_spindown(state.reference, state.solution.amplitudes, time_span, DECAY_CONSTANT_RATIO;
                                  initial_background=state.background, initial_guess=state.solution, newton_tol=newton_tolerance)
    return result, result.ode_solution.stats.nf
end

# --- run metadata ---------------------------------------------------------------

git_output(arguments...) = try
    readchomp(`git -C $(dirname(BENCHMARK_DIRECTORY)) $(arguments)`)
catch
    "unknown"
end

"""Commit hash, whether `src/` or `Project.toml` differ from it, host, CPU model and thread counts."""
function run_metadata()
    cpu_model = try
        strip(split(first(filter(line -> startswith(line, "Model name"), readlines(`lscpu`))), ':'; limit=2)[2])
    catch
        "unknown"
    end
    modified = git_output("status", "--porcelain", "--", "src", "Project.toml")
    return (commit=git_output("rev-parse", "--short=10", "HEAD"), source_modified=!isempty(modified),
            hostname=gethostname(), cpu_model, julia_version=string(VERSION), julia_threads=Threads.nthreads(),
            blas_threads=BLAS.get_num_threads(), timestamp=Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))
end

include("quantities.jl")

end # module
