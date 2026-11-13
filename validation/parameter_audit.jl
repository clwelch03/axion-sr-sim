# Numerical-parameter audit: convergence of the quantities one spindown right-hand-side evaluation produces
# (mode frequencies, d ln A/dt, db/dt, dj/dt) in one discretisation parameter at a time, around the production
# settings of scripts/evolution_spindown.jl, at the benchmark harness's states (benchmarks/common.jl) plus a
# takeover state. Results are written as CSV rows; `report` turns them into the tables of docs/parameter-audit.md.
#
#   julia --project=. -t 8 validation/parameter_audit.jl solve <study> <output-dir> [state,state,...]
#   julia --project=. -t 1 validation/parameter_audit.jl time <study> <output-dir> [state,...]   # reuses solved states
#   julia --project=. validation/parameter_audit.jl tolerance <output-dir> [state,...]
#   julia --project=. validation/parameter_audit.jl report <output-dir>
#
# Studies: see `STUDIES`. How a state is carried to a new discretisation: when the unknown layout is unchanged
# (quadrature changes only), the stored harness vector is the Newton guess; otherwise the state is rebuilt from
# scratch as benchmarks/build_states.jl does (cold solve for the linear states, cold solve at 1e-3 of the amplitudes
# and 30 geometric amplitude steps for the others). Every state is then solved to a scaled residual of 1e-12
# (floor_tol 1e-9 for meshes whose roundoff floor is higher), and the rates are evaluated with radial solves at the
# same tolerance, so that Newton error is far below the discretisation changes being measured.
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "..", "benchmarks", "common.jl"))
using .BenchmarkHarness
using AxionCloud
const AC = AxionCloud
const H = BenchmarkHarness

"""Discretisation and truncation settings varied by the audit; defaults are the production values."""
Base.@kwdef struct Discretisation
    n_phase::Int = 16
    n_polar::Int = 12
    step::Float64 = 0.4
    stretch_factor::Float64 = 1.03
    containment::Float64 = 1e-4
    horizon_offset_ratio::Float64 = 1e-6
    max_order::Int = 3
    n_angular::Int = 1
end

const PRODUCTION = Discretisation()
const FINER = Discretisation(n_phase=24, n_polar=16, step=0.3, stretch_factor=1.015, containment=1e-6,
                             horizon_offset_ratio=1e-8)

const STUDIES = Dict(
    "n_phase" => [Discretisation(n_phase=n) for n in (8, 12, 16, 24, 32)],
    "n_polar" => [Discretisation(n_polar=n) for n in (8, 10, 12, 16, 24)],
    "step" => [Discretisation(step=s) for s in (0.8, 0.6, 0.4, 0.3, 0.2)],
    "stretch" => [Discretisation(stretch_factor=f) for f in (1.06, 1.03, 1.015)],
    "outer" => [Discretisation(containment=c) for c in (1e-3, 1e-4, 1e-5, 1e-6, 1e-8)],
    "inner" => [Discretisation(horizon_offset_ratio=r) for r in (1e-4, 1e-5, 1e-6, 1e-8)],
    "truncation" => [Discretisation(), Discretisation(max_order=5), Discretisation(n_angular=2)],
    "combined" => [PRODUCTION, FINER],
)

const SOLVE_TOLERANCE = 1e-12
const FLOOR_TOLERANCE = 1e-9
# The α = 0.187 mesh has a residual floor between 2e-10 and 1e-9 that depends on the mesh and moves above 1e-9 after a
# perturbation; its solves accept a floor up to 3e-9.
floor_tolerance(spec) = spec.name == "large_mesh" ? 3e-9 : FLOOR_TOLERANCE

"""The takeover state: t = 2.94e10 of the restart run (|322⟩ rising fast while |211⟩ falls), not in the harness."""
const TAKEOVER = H.StateSpec("takeover", "|322> takeover, t = 2.9408e10 of log 61110567's solution CSV", 0.3, 0.9,
                             [2.551522594444469e-03, 1.950591472117652e-06], 9.840819661230717e-01, 8.429987750916837e-01,
                             :continuation, 1e-12, 1e-10, false)

"""
Two strongly nonlinear states of the f_a/M_Pl = 3e-3 runs (amplitudes about 7 times the plateau's), not in the harness:
the |211⟩ peak (t = 6.6667e8 of results/spindown_a0.9_alpha0.3_fa3e-3_t1e10_solution.csv) and the peak with |322⟩
grown (t = 6.75e8 of results/spindown_a0.9_alpha0.3_fa3e-3_solution.csv). Their rates use f_a/M_Pl = 3e-3.
"""
const STRONG_STATES = [
    H.StateSpec("strong_peak", "|211> peak of the f_a/M_Pl = 3e-3 run", 0.3, 0.9, [2.06565792e-2, 1.89863626e-6],
                0.9972886146555608, 0.8908392725156438, :continuation, 1e-12, 1e-10, false),
    H.StateSpec("strong_takeover", "both modes grown, f_a/M_Pl = 3e-3 run", 0.3, 0.9, [1.55687785e-2, 3.15544685e-5],
                0.9975424292397707, 0.8888775116704426, :continuation, 1e-12, 1e-10, false),
]

const STATE_NAMES = ["linear", "plateau", "takeover", "nonlinear", "large_mesh"]   # default; strong_* only when listed
function spec_of(name)
    name == "takeover" && return TAKEOVER
    index = findfirst(spec -> spec.name == name, STRONG_STATES)
    return index === nothing ? H.spec_named(name) : STRONG_STATES[index]
end
decay_constant_ratio(spec) = startswith(spec.name, "strong") ? 3e-3 : H.DECAY_CONSTANT_RATIO

"""Reference problem (initial background) of a state at the given discretisation."""
function reference_problem(spec, settings::Discretisation)
    background = KerrBackground(spec.initial_spin)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; max_order=settings.max_order, n_angular=settings.n_angular)
    mesh = RadialMesh(background, spec.mass_coupling, mode_set; step=settings.step, stretch_factor=settings.stretch_factor,
                      containment=settings.containment, horizon_offset_ratio=settings.horizon_offset_ratio)
    return RadialProblem(background, CosinePotential(), spec.mass_coupling, mode_set; mesh, n_polar=settings.n_polar,
                         n_phase=settings.n_phase)
end

same_layout(settings::Discretisation) =
    (settings.step, settings.stretch_factor, settings.containment, settings.horizon_offset_ratio, settings.max_order,
     settings.n_angular) == (PRODUCTION.step, PRODUCTION.stretch_factor, PRODUCTION.containment, PRODUCTION.horizon_offset_ratio,
                             PRODUCTION.max_order, PRODUCTION.n_angular)

"""Solve that throws on failure (caught by the caller, which records the failure)."""
function checked_solve(problem, amplitudes; guess=nothing, tol=SOLVE_TOLERANCE, floor_tol=FLOOR_TOLERANCE)
    solution = solve_radial(problem, amplitudes; guess, tol, floor_tol=max(tol, floor_tol))
    solution.converged || error("solve did not converge at amplitudes $amplitudes (|F| = $(solution.residual_norm))")
    return solution
end

"""
The converged state on `problem`: from the stored harness vector when there is one and the layout is unchanged, otherwise rebuilt from
scratch (cold solve, or cold solve at 1e-3 of the amplitudes plus 30 amplitude steps at tolerance 1e-10), then
solved to `SOLVE_TOLERANCE`. Returns the solution and the route taken.
"""
function converged_state(spec, settings::Discretisation, problem)
    if same_layout(settings) && isfile(H.state_path(spec))
        stored = collect(reinterpret(Float64, read(H.state_path(spec))))
        guess = H.solution_from_unknowns(problem, spec.amplitudes, stored, maximum(spec.amplitudes))
        return checked_solve(problem, spec.amplitudes; guess, floor_tol=floor_tolerance(spec)), "stored vector"
    end
    guess = nothing
    factors = spec.construction == :cold ? [1.0] : exp.(range(log(1e-3), 0.0; length=30))
    for factor in factors
        guess = checked_solve(problem, spec.amplitudes .* factor; guess, tol=max(1e-10, spec.solve_tolerance),
                              floor_tol=floor_tolerance(spec))
    end
    return checked_solve(problem, spec.amplitudes; guess, floor_tol=floor_tolerance(spec)), spec.construction == :cold ? "cold solve" : "continuation"
end

"""Background state of a spec and the problem there (current-mass units, mesh of the reference problem)."""
state_problem(spec, reference) = H.problem_at(reference, H.background_state(spec))

"""Evolution cache with production settings (action balance, no metric work) at the given Newton tolerances."""
audit_cache(spec, reference; tol, floor_tol) =
    EvolutionCache(reference; decay_constant_ratio=decay_constant_ratio(spec), newton_tol=tol, newton_floor_tol=floor_tol)

"""
Integrator-state rates `[d ln A₂₁₁/dt, d ln A₃₂₂/dt, db/dt, dj/dt]` from one `_spindown_rhs!` call at
`amplitudes` and `background`, warm-started from `solution` (no tangent prediction). Errors if the solve fails.
"""
function rhs_rates(cache, solution, amplitudes, background)
    state_vector = vcat(log.(amplitudes), background.mass, background.angular_momentum)
    cache.accepted_solution, cache.accepted_tangents, cache.failed_attempt = solution, nothing, -1
    rates = zeros(length(amplitudes) + 2)
    AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    cache.failed_attempt == -1 || error("right-hand-side solve failed")
    empty!(cache.diagnostics)
    return rates, cache.last_solution
end

"""The perturbed evaluation timed by the audit: the harness's fixed perturbation (A₃₂₂ × 1.01, b and j shifted)."""
perturbed_point(spec, solution) =
    (solution.amplitudes .* H.PERTURBATION_AMPLITUDE_FACTORS,
     BackgroundState(H.background_state(spec).mass * H.PERTURBATION_BACKGROUND_FACTORS[1],
                     H.background_state(spec).angular_momentum * H.PERTURBATION_BACKGROUND_FACTORS[2]))

"""Production Newton tolerances for a state (the α = 0.187 mesh needs the opt-in floor tolerance)."""
production_tolerances(spec) = (tol=1e-10, floor_tol=spec.name == "large_mesh" ? floor_tolerance(spec) : 1e-10)

"""Minimum wall time over `repeats` perturbed right-hand-side calls at production tolerances (after one warm-up)."""
function time_rhs(spec, reference, solution; repeats=3)
    cache = audit_cache(spec, reference; production_tolerances(spec)...)
    amplitudes, background = perturbed_point(spec, solution)
    _, warm = rhs_rates(cache, solution, amplitudes, background)
    best = minimum(@elapsed(rhs_rates(cache, solution, amplitudes, background)) for _ in 1:repeats)
    return best, warm.iterations
end

settings_label(settings::Discretisation) =
    @sprintf("np%d_npol%d_h%g_sf%g_c%g_hr%g_mo%d_na%d", settings.n_phase, settings.n_polar, settings.step, settings.stretch_factor,
             settings.containment, settings.horizon_offset_ratio, settings.max_order, settings.n_angular)

state_file(directory, study, settings, name) = joinpath(directory, "states", "$(study)_$(settings_label(settings))_$(name).jls")

const QUANTITY_NAMES = ["Re_omega_211", "Re_omega_322", "Im_omega_211", "Im_omega_322", "dlnA211_dt", "dlnA322_dt", "db_dt", "dj_dt"]

function append_row(path, fields)
    new_file = !isfile(path)
    open(path, "a") do io
        new_file && println(io, "study,label,state,kind,threads,nodes,unknowns,residual,route,seconds,newton," * join(QUANTITY_NAMES, ","))
        println(io, join(fields, ","))
    end
end

format(value::Real) = @sprintf("%.15e", value)
format(value) = string(value)

"""Solve each state at each setting of the study, evaluate the quantities, and time the perturbed call at this thread count."""
function run_solve(study, directory, names)
    mkpath(joinpath(directory, "states"))
    csv = joinpath(directory, "$study.csv")
    for settings in STUDIES[study], name in names
        spec = spec_of(name)
        label = settings_label(settings)
        try
            reference = reference_problem(spec, settings)
            problem = state_problem(spec, reference)
            build_seconds = @elapsed (solution, route) = converged_state(spec, settings, problem)
            serialize(state_file(directory, study, settings, name), solution.unknowns)
            frequencies = mode_frequencies(solution)
            cache = audit_cache(spec, reference; tol=SOLVE_TOLERANCE, floor_tol=floor_tolerance(spec))
            rates, _ = rhs_rates(cache, solution, solution.amplitudes, H.background_state(spec))
            quantities = vcat(real.(frequencies), imag.(frequencies), rates)
            append_row(csv, [study, label, name, "quantities", Threads.nthreads(), length(reference.mesh.rstar),
                             length(solution.unknowns), format(solution.residual_norm), route, format(build_seconds), solution.iterations,
                             format.(quantities)...])
            seconds, iterations = time_rhs(spec, reference, solution)
            append_row(csv, [study, label, name, "timing", Threads.nthreads(), length(reference.mesh.rstar), length(solution.unknowns),
                             "", "", format(seconds), iterations, fill("", length(QUANTITY_NAMES))...])
            @printf("%-10s %-50s %-10s %5d nodes  %-13s build %6.1f s  rhs %.3f s (%d it)  Im ω = (%.10e, %.10e)  rates = %s\n",
                    study, label, name, length(reference.mesh.rstar), route, build_seconds, seconds, iterations,
                    imag(frequencies[1]), imag(frequencies[2]), join((@sprintf("%.10e", r) for r in rates), ", "))
            name == "plateau" && run_segment(csv, study, label, spec, reference, solution)
        catch failure
            message = replace(sprint(showerror, failure), r"[\n,]" => " ")
            append_row(csv, [study, label, name, "failed", Threads.nthreads(), "", "", "", message[1:min(end, 200)], "", "",
                             fill("", length(QUANTITY_NAMES))...])
            println("$study $label $name FAILED: $message")
        end
        flush(stdout)
    end
end

"""
The harness's evolution segment 1.5e9 → 1.9e9 from the plateau state at this discretisation (Tsit5, library default
tolerances, Newton tolerance 1e-10). Row `segment`: final [ln A₂₁₁, ln A₃₂₂, b, j] in the last four quantity
columns (the first four hold the right-hand-side call count, accepted steps, and zeros), wall time in `seconds`.
"""
function run_segment(csv, study, label, spec, reference, solution)
    seconds = @elapsed result = evolve_with_spindown(reference, solution.amplitudes, (1.5e9, 1.9e9), H.DECAY_CONSTANT_RATIO;
                                                     initial_background=H.background_state(spec), initial_guess=solution,
                                                     newton_tol=1e-10)
    final_state = result.ode_solution.u[end]
    stats = result.ode_solution.stats
    append_row(csv, [study, label, spec.name, "segment", Threads.nthreads(), length(reference.mesh.rstar), length(solution.unknowns),
                     "", "", format(seconds), "", stats.nf, stats.naccept, 0, 0, format.(final_state)...])
    @printf("%-10s %-50s segment %.2f s, %d rhs calls, final state %s\n", study, label, seconds, stats.nf,
            join((@sprintf("%.12e", x) for x in final_state), ", "))
end

"""Time the perturbed call for states solved earlier by `run_solve` (at this run's thread count)."""
function run_time(study, directory, names)
    csv = joinpath(directory, "$study.csv")
    for settings in STUDIES[study], name in names
        path = state_file(directory, study, settings, name)
        isfile(path) || continue
        spec = spec_of(name)
        reference = reference_problem(spec, settings)
        problem = state_problem(spec, reference)
        unknowns = deserialize(path)
        solution = H.solution_from_unknowns(problem, spec.amplitudes, unknowns, maximum(spec.amplitudes))
        seconds, iterations = time_rhs(spec, reference, solution)
        append_row(csv, [study, settings_label(settings), name, "timing", Threads.nthreads(), length(reference.mesh.rstar),
                         length(unknowns), "", "", format(seconds), iterations, fill("", length(QUANTITY_NAMES))...])
        @printf("%-10s %-50s %-10s rhs %.3f s (%d it) at %d threads\n", study, settings_label(settings), name, seconds, iterations,
                Threads.nthreads())
        flush(stdout)
    end
end

"""
Newton tolerance study at production discretisation: the perturbed evaluation solved from the converged state at
`tol` (radial solves of the rates at the same `tol`), against the same evaluation at 1e-12.
"""
function run_tolerance(directory, names)
    mkpath(directory)
    csv = joinpath(directory, "tolerance.csv")
    for name in names
        spec = spec_of(name)
        reference = reference_problem(spec, PRODUCTION)
        problem = state_problem(spec, reference)
        solution, route = converged_state(spec, PRODUCTION, problem)
        amplitudes, background = perturbed_point(spec, solution)
        for tol in (1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12)
            try
                cache = audit_cache(spec, reference; tol, floor_tol=max(tol, floor_tolerance(spec)))
                rates, evaluated = rhs_rates(cache, solution, amplitudes, background)
                seconds = minimum(@elapsed(rhs_rates(cache, solution, amplitudes, background)) for _ in 1:3)
                frequencies = mode_frequencies(evaluated)
                append_row(csv, ["tolerance", @sprintf("tol%g", tol), name, "quantities", Threads.nthreads(), length(reference.mesh.rstar),
                                 length(evaluated.unknowns), format(evaluated.residual_norm), route, format(seconds), evaluated.iterations,
                                 format.(vcat(real.(frequencies), imag.(frequencies), rates))...])
                @printf("tolerance %g %-10s |F| = %.2e  %d it  %.3f s  rates = %s\n", tol, name, evaluated.residual_norm,
                        evaluated.iterations, seconds, join((@sprintf("%.12e", r) for r in rates), ", "))
            catch failure
                println("tolerance $tol $name FAILED: ", sprint(showerror, failure))
            end
            flush(stdout)
        end
    end
end

function main(arguments)
    command = arguments[1]
    if command == "solve" || command == "time"
        study, directory = arguments[2], arguments[3]
        names = length(arguments) ≥ 4 ? split(arguments[4], ",") : STATE_NAMES
        command == "solve" ? run_solve(study, directory, names) : run_time(study, directory, names)
    elseif command == "tolerance"
        run_tolerance(arguments[2], length(arguments) ≥ 3 ? split(arguments[3], ",") : STATE_NAMES)
    else
        error("unknown command $command")
    end
end

main(ARGS)
