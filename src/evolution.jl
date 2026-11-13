"""
Adiabatic time evolution of the cloud amplitudes.

At every instant the cloud is the stationary multiphase solution at the current
amplitudes; a balance between the slowly changing conserved quantities (the
wave actions by default, or energy and angular momentum) and their boundary
fluxes then fixes the amplitude rates. Two levels are provided:

* `evolve_fixed_background`: the black hole does not respond (toy model),
  `[∂_A N] A' = −L_N`.
* `evolve_with_spindown`: the black-hole mass and spin absorb the horizon
  fluxes, and the balance includes the background derivatives of the charges
  and the rate of metric work on the cloud.

Two equivalent balances are available (`balance` keyword):

* `:action` (default since 2026-11-09): one wave-action balance per independent phase,
  `[∂_A N] A' = −L_N − (∂_b N b' + ∂_j N j')`, valid for any number of modes. The
  actions are adiabatic invariants, so no metric-work term appears, and their
  fluxes contain no time derivative. Energy and angular momentum follow from
  `J = Σ m_i N_i` and `dE = Σ ω_i dN_i + w`. On a fixed background the two
  balances agree up to the `Im ω̃ Re(R*Q)` term that `outward_flux` omits from
  the energy flux at a boundary where the cloud is not contained
  (`validation/action_balance.jl`).
* `:energy` (the outline's formulation, the default before 2026-11-09): energy and
  angular momentum, `[∂_A E; ∂_A J] A' = −L + (w, 0) − (∂_b Q b' + ∂_j Q j')`. Two
  conservation laws, so exactly two bound modes. Pass `balance=:energy` to
  reproduce results computed before the default changed.

Reference units are the *initial* black-hole mass `M₀` (and `f_a`): with
`b = M/M₀`, `j = J_BH/M₀²`, `τ = t/M₀`, a radial solution in units of the
current mass carries `E_ref = b E`, `J_ref = b² J`, energy losses per unit `τ`
equal to the current-unit fluxes, and angular-momentum losses `b` times them;
actions scale like `J` (`N_ref = b² N`, losses `b` times the current-unit ones).
The action balance handles any number of bound modes; the energy balance (two
conservation laws, two amplitude equations) requires exactly two.

Cloud charges are measured in units of `f_a² M₀` while the black hole's are in
units of `M₀` (with `G = 1`), so the horizon fluxes change the hole's mass and
spin only after multiplication by `(f_a / M_Pl)²`, where `M_Pl = G^{-1/2}` is the
(non-reduced) Planck mass. This `decay_constant_ratio = f_a / M_Pl` is a
required input of the spindown driver.

Both drivers integrate the *logarithms* of the amplitudes by default
(`log_amplitudes=true`): in the linear regime `d ln A/dt = Im ω` is constant, so
the step size is limited by how fast the rates change rather than by the
exponential growth itself, and the error tolerances become relative accuracies
of the amplitudes. The first step is set from the initial e-folding times
(`initial_step_fraction`) rather than by the integrator's heuristic, which the
`1e-13`-level noise of the Newton solve otherwise drives down to `dt ≈ 1`.
Use `amplitude_trajectory` to read the amplitudes back from a result.
"""

"""
    EvolutionDiagnostics

Snapshot recorded at every right-hand-side evaluation: time, amplitudes,
complex mode frequencies, reference-unit charges and losses, and (with
spindown) the background state and metric work; also the occupations
`N_i` (`phase_actions`) and their losses `L_{N_i}` (`action_losses`), in
reference units.
"""
struct EvolutionDiagnostics{T<:Real}
    time::T
    amplitudes::Vector{T}
    mode_frequencies::Vector{Complex{T}}
    charges::Charges{T}
    losses::Charges{T}
    black_hole_mass::T
    black_hole_spin::T
    metric_work::T
    actions::Vector{T}
    action_losses::Vector{T}
end

"""
    RadialSolveCounts

Radial solves of one evolution by purpose: `stage` (right-hand-side evaluations
requested by the integrator, the entries of the diagnostics log), `jacobian`
(perturbed solves of the finite-difference Jacobian), `accepted_state` (the
initial solve, and the solves at accepted states that a new Jacobian needs when
the integrator's last stage is not there), and `failed` (solves of any purpose
that did not converge, included in the other counts). `jacobians_computed` and
`jacobians_reused` count Jacobian requests by outcome.
"""
mutable struct RadialSolveCounts
    stage::Int
    jacobian::Int
    accepted_state::Int
    failed::Int
    jacobians_computed::Int
    jacobians_reused::Int
end

RadialSolveCounts() = RadialSolveCounts(0, 0, 0, 0, 0, 0)

"""Total number of radial solves recorded in `counts`."""
total_solves(counts::RadialSolveCounts) = counts.stage + counts.jacobian + counts.accepted_state

"""
    EvolutionResult

Output of a driver: the `OrdinaryDiffEq` solution object (state `[A₁, …, A_K]` or
`[A₁, …, A_K, b, j]`, with `ln A` in place of `A` when `log_amplitudes` is true),
the diagnostics log (one entry per integrator stage evaluation), the final
radial solution, and the radial solves by purpose (`RadialSolveCounts`). Read
the amplitudes with `amplitude_trajectory` rather than from `ode_solution.u`
directly.
"""
struct EvolutionResult{S, T<:Real, R<:RadialSolution}
    ode_solution::S
    diagnostics::Vector{EvolutionDiagnostics{T}}
    final_radial_solution::R
    log_amplitudes::Bool
    solve_counts::RadialSolveCounts
end

"""
    AmplitudeTrajectory

Amplitudes at the integrator's output times, whatever the state convention:
`times`, `amplitudes[time_index, mode]`, and the black-hole mass `b` and
angular momentum `j` (constant `1` and `a` on a fixed background).
"""
struct AmplitudeTrajectory{T<:Real}
    times::Vector{T}
    amplitudes::Matrix{T}
    black_hole_mass::Vector{T}
    black_hole_angular_momentum::Vector{T}
end

"""
    amplitude_trajectory(result)

Amplitudes (and background state) at the saved times of an evolution, undoing
the logarithmic state if it was used.
"""
function amplitude_trajectory(result::EvolutionResult)
    states = result.ode_solution.u
    mode_count = n_modes(result.final_radial_solution.problem)
    transform = result.log_amplitudes ? exp : identity
    amplitudes = permutedims(reduce(hcat, [transform.(state[1:mode_count]) for state in states]))
    with_background = length(first(states)) == mode_count + 2
    spin = result.final_radial_solution.problem.background.a
    mass = with_background ? [state[mode_count + 1] for state in states] : ones(length(states))
    angular_momentum = with_background ? [state[mode_count + 2] for state in states] : fill(spin, length(states))
    AmplitudeTrajectory(collect(result.ode_solution.t), amplitudes, mass, angular_momentum)
end

"""
    SolutionTangents

First-order response ("tangents") of a converged `solution`'s unknown vector, in units of its `scale` with the
frequencies unscaled: to each prescribed amplitude (columns of `wrt_amplitudes`, from `amplitude_response`), and at
fixed amplitudes to the black-hole mass `b` and angular momentum `j` (`background_tangents`), at the
background `(black_hole_mass, black_hole_angular_momentum)` of the solution. Built by `solution_tangents`; used by
`predicted_change` to start a nearby solve within second order of its solution.
"""
struct SolutionTangents{T<:Real, S<:RadialSolution}
    solution::S
    black_hole_mass::T
    black_hole_angular_momentum::T
    wrt_amplitudes::Matrix{T}
    wrt_mass::Vector{T}
    wrt_angular_momentum::Vector{T}
end

"""Mutable scratch shared by the right-hand-side evaluations of one evolution."""
mutable struct EvolutionCache{P<:RadialProblem, T<:Real}
    template::P
    initial_mass_coupling::T
    coupling_to_hole::T   # (f_a / M_Pl)²: converts cloud fluxes in f_a² units into black-hole mass rates
    include_metric_work::Bool
    record_metric_work::Bool   # also compute E, J and the metric work under the action balance (diagnostic only)
    last_solution::Union{Nothing, RadialSolution}
    diagnostics::Vector{EvolutionDiagnostics{T}}
    newton_tol::T
    verbose::Bool
    log_amplitudes::Bool
    balance::Symbol       # :energy or :action
    accepted_solution::Union{Nothing, RadialSolution}   # radial solution at the last accepted step (spindown warm start)
    solve_failed::Bool    # a radial solve in the current step attempt did not converge (spindown rejects the step)
    tangent_predictor::Bool   # start solves from a first-order prediction (`predicted_change`, `_amplitude_prediction`)
    last_rates::Vector{T}       # integrator-state rates at `last_solution` (spindown)
    accepted_rates::Vector{T}   # integrator-state rates at `accepted_solution` (spindown)
    last_tangents::Union{Nothing, SolutionTangents}       # tangents of `last_solution` (spindown with tangent_predictor)
    accepted_tangents::Union{Nothing, SolutionTangents}   # tangents of `accepted_solution`
    solve_counts::RadialSolveCounts
end

"""
    EvolutionCache(problem; decay_constant_ratio=0.0, include_metric_work=false, record_metric_work=false,
                   newton_tol=1e-10, verbose=false, log_amplitudes=true, balance=:action, last_solution=nothing,
                   tangent_predictor=true)

Fresh cache for an evolution of `problem` (also handy for calling
`balance_terms` directly). The energy, the angular momentum and the metric work
are computed only when the energy balance, `include_metric_work` or
`record_metric_work` asks for them (`_computes_charges`).
"""
function EvolutionCache(problem::RadialProblem; decay_constant_ratio=0.0, include_metric_work=false, record_metric_work=false,
                        newton_tol=1e-10, verbose=false, log_amplitudes=true, balance=:action, last_solution=nothing,
                        tangent_predictor=true)
    _check_balance(problem, balance)
    EvolutionCache(problem, problem.mass_coupling, float(decay_constant_ratio)^2, include_metric_work, record_metric_work, last_solution,
                   EvolutionDiagnostics{Float64}[], float(newton_tol), verbose, log_amplitudes, balance,
                   last_solution, false, tangent_predictor, Float64[], Float64[], nothing, nothing, RadialSolveCounts())
end

"""Whether the balance terms include the energy, angular momentum and metric work: the action balance does not need them."""
_computes_charges(cache::EvolutionCache) = cache.balance == :energy || cache.include_metric_work || cache.record_metric_work

"""Print one progress line per right-hand-side evaluation when the cache is verbose."""
function _report_step(cache::EvolutionCache, time, amplitudes, solution::RadialSolution, state_text::String)
    cache.verbose || return nothing
    println("t = ", time, "  A = ", amplitudes, "  Im ω = ", growth_rates(solution), "  newton = ", solution.iterations, state_text)
    flush(stdout)
    return nothing
end

"""Reject unknown balances, and the energy balance for anything but two bound modes."""
function _check_balance(problem::RadialProblem, balance::Symbol)
    balance in (:energy, :action) || throw(ArgumentError("balance must be :energy or :action; got :$balance"))
    balance == :energy && n_modes(problem) != 2 &&
        throw(ArgumentError("the energy balance needs exactly two bound modes (energy and angular momentum give two " *
                            "equations); got $(n_modes(problem)). Use balance=:action for any number of modes"))
    return nothing
end

"""Amplitudes from the integrator state's first `n_modes` entries (`ln A` or `A`)."""
_state_amplitudes(state, n_modes::Int, log_amplitudes::Bool) =
    log_amplitudes ? exp.(state[1:n_modes]) : state[1:n_modes]

"""Integrator-state rates of the amplitude entries from the amplitude rates."""
_state_rates(amplitude_rates, amplitudes, log_amplitudes::Bool) = log_amplitudes ? amplitude_rates ./ amplitudes : amplitude_rates

"""
Initial step `initial_step_fraction / max_i |A_i'/A_i|` from the rates at the
start, or `nothing` (the integrator's own heuristic) if every rate vanishes or
`initial_step_fraction` is `nothing`.
"""
function _initial_step(amplitudes, amplitude_rates, initial_step_fraction)
    initial_step_fraction === nothing && return nothing
    fastest = maximum(abs.(amplitude_rates ./ amplitudes))
    return fastest > 0 ? initial_step_fraction / fastest : nothing
end

"""Solver keywords with `dt` from `_initial_step` unless the caller supplied one."""
function _with_initial_step(kwargs, step)
    (haskey(kwargs, :dt) || step === nothing) && return kwargs
    return merge(NamedTuple(kwargs), (dt=step,))
end

"""
First-order change of `previous`'s unknowns for moving to `amplitudes` on the same discrete problem
(`amplitude_response`, one back-substitution per mode), or `nothing` when the tangent predictor is off or
`previous` is not a solution of `problem`.
"""
function _amplitude_prediction(cache::EvolutionCache, problem::RadialProblem, previous, amplitudes)
    (cache.tangent_predictor && previous isa RadialSolution && _same_discretisation(problem, previous.problem)) || return nothing
    return amplitude_response(previous) * (amplitudes .- previous.amplitudes)
end

"""Radial solve warm-started from the cache's previous solution, moved by `_amplitude_prediction`."""
function _solve_cached!(cache::EvolutionCache, problem::RadialProblem, amplitudes)
    previous = cache.last_solution
    solution = solve_radial(problem, amplitudes; guess=previous,
                            guess_change=_amplitude_prediction(cache, problem, previous, amplitudes), tol=cache.newton_tol)
    cache.last_solution = solution
    return solution
end

"""Print one line that is not a stage evaluation (Jacobian, accepted-state solve) when the cache is verbose."""
function _report_line(cache::EvolutionCache, text::String)
    cache.verbose || return nothing
    println(text)
    flush(stdout)
    return nothing
end

"""
Whether `solution` was solved at the integrator state `state`: same amplitudes
and, when the state carries `b` and `j`, the same spin `j/b²` and mass coupling
`α₀ b`, all to `rtol = 1e-12`.
"""
function _solution_at_state(solution, state, cache::EvolutionCache)
    solution === nothing && return false
    mode_count = n_modes(cache.template)
    isapprox(solution.amplitudes, _state_amplitudes(state, mode_count, cache.log_amplitudes); rtol=1e-12) || return false
    length(state) == mode_count && return true
    mass, angular_momentum = state[mode_count + 1], state[mode_count + 2]
    return isapprox(solution.problem.background.a, angular_momentum / mass^2; rtol=1e-12) &&
           isapprox(solution.problem.mass_coupling, cache.initial_mass_coupling * mass; rtol=1e-12)
end

# --- finite-difference Jacobian for stiff integrators ----------------------------

"""
    FiniteDifferenceJacobian(; difference=:forward, amplitude_step=3e-4, mass_step=1e-6,
                             angular_momentum_step=1e-6, step_reductions=2, max_reuse_steps=0)

Settings of the Jacobian of the right-hand side with respect to the integrator
state, which the drivers supply to implicit and Rosenbrock integrators
(`algorithm`) in place of the integrator's own differencing. Column `k` is a
`:forward` or `:central` difference quotient with an absolute step in state
entry `k`: `amplitude_step` in `ln A` (a step relative to `A` without
`log_amplitudes`), `mass_step` in `b`, `angular_momentum_step` in `j`. Every
perturbed radial solve starts from the solution at the base state, and the
perturbed solves run as concurrent tasks (`Threads.@spawn`). A perturbed solve
that fails is retried with its step divided by 4, at most `step_reductions`
times. `max_reuse_steps` is the reuse policy of `_reuse_stored_jacobian`: the
default `0` computes a new Jacobian at every accepted state.

The default steps come from step scans at a = 0.9, α = 0.3, f_a/M_Pl = 0.047
in the linear-growth, |211⟩-plateau and |322⟩-takeover regimes: there forward
differences agree with central ones to 2e-4 (amplitude columns) and 3e-5 (`b`
and `j` columns), the noise of the rates contributes less than 1e-3, and the
perturbed solves stay clear of the outliers seen with `b` steps of 1e-7 (rates
from sensitivities factorised away from the converged point, see `solve_radial`).
"""
struct FiniteDifferenceJacobian
    difference::Symbol
    amplitude_step::Float64
    mass_step::Float64
    angular_momentum_step::Float64
    step_reductions::Int
    max_reuse_steps::Int
end

function FiniteDifferenceJacobian(; difference::Symbol=:forward, amplitude_step=3e-4, mass_step=1e-6,
                                  angular_momentum_step=1e-6, step_reductions=2, max_reuse_steps=0)
    difference in (:forward, :central) || throw(ArgumentError("difference must be :forward or :central; got :$difference"))
    FiniteDifferenceJacobian(difference, amplitude_step, mass_step, angular_momentum_step, step_reductions, max_reuse_steps)
end

"""
Jacobian kept between the integrator's requests: the matrix, the state and
time at which it was computed (`computed_state` empty before the first), the
state of the previous request, and its age, the number of distinct states
requested since it was computed.
"""
mutable struct JacobianStore
    matrix::Matrix{Float64}
    computed_state::Vector{Float64}
    computed_time::Float64
    last_request_state::Vector{Float64}
    age::Int
end

JacobianStore(state_length::Int) = JacobianStore(zeros(state_length, state_length), Float64[], NaN, Float64[], 0)

"""
Whether the stored Jacobian serves a request at `state`; updates the request
bookkeeping. The integrators request the Jacobian at the start of every step
attempt, at the last accepted state, so a request at the same state as the
previous one means the previous attempt was rejected. A request at the state
where the stored Jacobian was computed always reuses it (it is exact there).
Otherwise the stored Jacobian is reused while at most `max_reuse_steps` new
states have been requested since it was computed, except when the request
repeats the previous one: a step rejected with a stale Jacobian is retried with
a new one.
"""
function _reuse_stored_jacobian(store::JacobianStore, settings::FiniteDifferenceJacobian, state::AbstractVector)
    repeated = state == store.last_request_state
    if !repeated
        store.age += 1
        store.last_request_state = copy(state)
    end
    isempty(store.computed_state) && return false
    state == store.computed_state && return true
    return !repeated && store.age <= settings.max_reuse_steps
end

"""Absolute difference steps of the state entries (`ln A` or `A`, then `b` and `j` when present)."""
function _difference_steps(settings::FiniteDifferenceJacobian, state::AbstractVector, mode_count::Int, log_amplitudes::Bool)
    amplitude_steps = log_amplitudes ? fill(settings.amplitude_step, mode_count) :
                                       settings.amplitude_step .* abs.(state[1:mode_count])
    length(state) == mode_count && return amplitude_steps
    return vcat(amplitude_steps, settings.mass_step, settings.angular_momentum_step)
end

"""
Evaluation at the base state plus `step` in entry `component`, through
`evaluate(perturbation)` (a named tuple whose `rates` is `nothing` when the
radial solve failed). A failed solve is retried with the step divided by 4, at
most `step_reductions` times. Returns the last evaluation, the step it used and
the number of solves attempted.
"""
function _perturbed_evaluation(evaluate, state_length::Int, component::Int, step, step_reductions::Int)
    evaluation, trial_step = nothing, step
    for attempt in 0:step_reductions
        trial_step = step / 4^attempt
        perturbation = zeros(state_length)
        perturbation[component] = trial_step
        evaluation = evaluate(perturbation)
        evaluation.rates === nothing || return (; evaluation, step=trial_step, attempts=attempt + 1)
    end
    return (; evaluation, step=trial_step, attempts=step_reductions + 1)
end

"""
Finite-difference Jacobian of the rates at `state`, with `base_rates` the rates
there and `evaluate` as in `_perturbed_evaluation`. The perturbed evaluations
(one per column for forward differences, two for central) run as concurrent
tasks; `evaluate` must therefore not mutate shared state. A central column
whose two sides ended with different steps uses the quotient over their
separation. Returns the matrix, or `nothing` if some column has a failed side,
and the perturbed results `[side, column]`.
"""
function _difference_jacobian(evaluate, state::AbstractVector, base_rates::AbstractVector, steps::AbstractVector,
                              settings::FiniteDifferenceJacobian)
    signs = settings.difference == :central ? (1, -1) : (1,)
    tasks = [Threads.@spawn _perturbed_evaluation(evaluate, length(state), component, sign * steps[component],
                                                  settings.step_reductions)
             for sign in signs, component in eachindex(state)]
    results = fetch.(tasks)
    any(result -> result.evaluation.rates === nothing, results) && return nothing, results
    matrix = zeros(length(base_rates), length(state))
    for component in eachindex(state)
        plus = results[1, component]
        lower_rates, lower_step = settings.difference == :central ?
            (results[2, component].evaluation.rates, results[2, component].step) : (base_rates, 0.0)
        matrix[:, component] = (plus.evaluation.rates .- lower_rates) ./ (plus.step - lower_step)
    end
    return matrix, results
end

"""
Fill `jacobian_matrix` for a request at `state` and record the solves. Reuses
the stored Jacobian when `_reuse_stored_jacobian` allows; otherwise `base()`
returns the radial solution and rates at `state` and `evaluate(base_solution,
perturbation)` the evaluation at `state + perturbation`. If a column cannot be
formed even with reduced steps, the stored Jacobian is used if there is one
(a W-type method such as `Rosenbrock23` tolerates it; rejecting the step would
not help, since the retry asks for the Jacobian at the same state); without a
stored Jacobian the evolution stops with an error.
"""
function _fill_jacobian!(jacobian_matrix, state, time, cache::EvolutionCache, settings::FiniteDifferenceJacobian,
                         store::JacobianStore, base, evaluate)
    counts = cache.solve_counts
    if _reuse_stored_jacobian(store, settings, state)
        counts.jacobians_reused += 1
        jacobian_matrix .= store.matrix
        _report_line(cache, "jacobian t = $time: reused (computed at t = $(store.computed_time), age $(store.age))")
        return nothing
    end
    base_solution, base_rates = base()
    steps = _difference_steps(settings, state, n_modes(cache.template), cache.log_amplitudes)
    matrix, results = _difference_jacobian(perturbation -> evaluate(base_solution, perturbation), state, base_rates,
                                           steps, settings)
    attempts = sum(result.attempts for result in results)
    counts.jacobian += attempts
    counts.failed += attempts - count(result -> result.evaluation.rates !== nothing, results)
    newton_iterations = [result.evaluation.solution.iterations for result in results]
    if matrix === nothing
        isempty(store.computed_state) &&
            error("finite-difference Jacobian at t = $time failed: perturbed radial solves did not converge with " *
                  "steps $steps reduced $(settings.step_reductions) times, and no earlier Jacobian is stored")
        jacobian_matrix .= store.matrix
        _report_line(cache, "jacobian t = $time: perturbed solves failed ($attempts solves); using the Jacobian " *
                            "computed at t = $(store.computed_time)")
        return nothing
    end
    counts.jacobians_computed += 1
    store.matrix, store.computed_state, store.computed_time, store.age = matrix, copy(state), time, 0
    jacobian_matrix .= matrix
    _report_line(cache, "jacobian t = $time: computed from $attempts solves, newton = $(vec(newton_iterations)), " *
                        "steps = $([result.step for result in results][1, :])")
    return nothing
end

"""Time gradient of the autonomous right-hand side (zero), so that stiff integrators do not difference in `t`."""
function _zero_time_gradient!(time_gradient, state, cache, time)
    time_gradient .= 0
    return nothing
end

"""ODE function with the Jacobian `jacobian!(matrix, state, cache, time)` and a zero time gradient."""
_ode_function(rhs!, jacobian!) = ODEFunction(rhs!; jac=jacobian!, tgrad=_zero_time_gradient!)

# --- fixed background -----------------------------------------------------------

"""
    amplitude_rates(solution; balance=:action)

`dA/dt` on a fixed background, from `[∂_A N] A' = −L_N` (`balance=:action`,
the default, any number of modes) or `[∂_A E; ∂_A J] A' = −[L_E; L_J]`
(`balance=:energy`, two modes).
"""
function amplitude_rates(solution::RadialSolution; balance::Symbol=:action)
    _check_balance(solution.problem, balance)
    balance == :action && return -(action_jacobian(solution) \ action_losses(solution))
    return -(charge_jacobian(solution) \ collect(boundary_losses(solution)))
end

"""Integrator-state rates on a fixed background from a radial solution at `amplitudes`."""
_fixed_background_rates(solution::RadialSolution, amplitudes, cache::EvolutionCache) =
    _state_rates(amplitude_rates(solution; cache.balance), amplitudes, cache.log_amplitudes)

"""
Radial solve on the fixed background warm-started from `guess`, and the rates at
the integrator state `state` (`nothing` if the solve did not converge). Mutates
nothing.
"""
function _fixed_background_evaluation(state, cache::EvolutionCache, guess)
    amplitudes = _state_amplitudes(state, n_modes(cache.template), cache.log_amplitudes)
    solution = solve_radial(cache.template, amplitudes; guess,
                            guess_change=_amplitude_prediction(cache, cache.template, guess, amplitudes), tol=cache.newton_tol)
    rates = solution.converged ? _fixed_background_rates(solution, amplitudes, cache) : nothing
    return (; rates, solution)
end

"""
Jacobian request on the fixed background. The base is the last stage solution
when it is at `state` (first-same-as-last integrators), otherwise one more solve.
"""
function _fixed_background_jacobian!(jacobian_matrix, state, cache::EvolutionCache, settings::FiniteDifferenceJacobian,
                                     store::JacobianStore, time)
    function base()
        if _solution_at_state(cache.last_solution, state, cache)
            amplitudes = _state_amplitudes(state, n_modes(cache.template), cache.log_amplitudes)
            return cache.last_solution, _fixed_background_rates(cache.last_solution, amplitudes, cache)
        end
        evaluation = _fixed_background_evaluation(state, cache, cache.last_solution)
        cache.solve_counts.accepted_state += 1
        evaluation.rates === nothing && error("radial solve at the Jacobian's base state t = $time did not converge")
        return evaluation.solution, evaluation.rates
    end
    _fill_jacobian!(jacobian_matrix, state, time, cache, settings, store, base,
                    (base_solution, perturbation) -> _fixed_background_evaluation(state .+ perturbation, cache, base_solution))
end

function _fixed_background_rhs!(rates, state, cache::EvolutionCache, time)
    amplitudes = _state_amplitudes(state, n_modes(cache.template), cache.log_amplitudes)
    solution = _solve_cached!(cache, cache.template, amplitudes)
    cache.solve_counts.stage += 1
    rates .= _fixed_background_rates(solution, amplitudes, cache)
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(amplitudes), mode_frequencies(solution),
                                                  cloud_charges(solution), boundary_losses(solution),
                                                  1.0, cache.template.background.a, 0.0,
                                                  phase_actions(solution), action_losses(solution)))
    _report_step(cache, time, amplitudes, solution, "")
    return nothing
end

"""
    evolve_fixed_background(problem, amplitudes, time_span; balance=:action, log_amplitudes=true,
                            initial_step_fraction=0.05, newton_tol=1e-10, verbose=false, tangent_predictor=true,
                            algorithm=Tsit5(), jacobian_settings=FiniteDifferenceJacobian(), reltol, abstol, kwargs...)

Integrate the amplitudes with `algorithm` (default Tsit5) on a frozen Kerr
background, with the action (default, any number of modes) or energy (two
modes) balance. Implicit and Rosenbrock algorithms get the finite-difference
Jacobian described by `jacobian_settings` (`FiniteDifferenceJacobian`).
With `log_amplitudes` (default) the state is `ln A` and the default tolerances
`reltol = abstol = 1e-7` bound the per-step relative error of the amplitudes
by about `1e-7 (1 + |ln A|)`; otherwise the state is `A` with the old defaults
`reltol = 1e-6`, `abstol = 1e-12`. The first step is `initial_step_fraction`
times the shortest initial e-folding time unless `dt` is passed
(`initial_step_fraction=nothing` leaves it to Tsit5's heuristic).
`verbose` prints one line per right-hand-side evaluation. Each radial solve is
warm-started from the previous one, moved to first order in the amplitude change
unless `tangent_predictor=false`. Other keyword arguments are passed to `solve`
(e.g. `saveat`, `dtmax`).
"""
function evolve_fixed_background(problem::RadialProblem, amplitudes::AbstractVector, time_span; balance::Symbol=:action,
                                 log_amplitudes=true, initial_step_fraction=0.05, newton_tol=1e-10, verbose=false,
                                 tangent_predictor=true, algorithm=Tsit5(),
                                 jacobian_settings::FiniteDifferenceJacobian=FiniteDifferenceJacobian(),
                                 reltol=(log_amplitudes ? 1e-7 : 1e-6), abstol=(log_amplitudes ? 1e-7 : 1e-12), kwargs...)
    cache = EvolutionCache(problem; newton_tol, verbose, log_amplitudes, balance, tangent_predictor)
    initial_amplitudes = collect(float.(amplitudes))
    initial_rates = amplitude_rates(_solve_cached!(cache, problem, initial_amplitudes); balance)
    cache.solve_counts.accepted_state += 1
    solver_options = _with_initial_step(kwargs, _initial_step(initial_amplitudes, initial_rates, initial_step_fraction))
    initial_state = log_amplitudes ? log.(initial_amplitudes) : initial_amplitudes
    store = JacobianStore(length(initial_state))
    jacobian! = (matrix, state, cache, time) -> _fixed_background_jacobian!(matrix, state, cache, jacobian_settings, store, time)
    ode = ODEProblem(_ode_function(_fixed_background_rhs!, jacobian!), initial_state, time_span, cache)
    ode_solution = solve(ode, algorithm; reltol, abstol, solver_options...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution, log_amplitudes, cache.solve_counts)
end

# --- evolving background --------------------------------------------------------

"""
    BackgroundState(mass, angular_momentum)

Black-hole mass `b = M/M₀` and angular momentum `j = J/M₀²` in reference units.
"""
struct BackgroundState{T<:Real}
    mass::T
    angular_momentum::T
end

dimensionless_spin(state::BackgroundState) = state.angular_momentum / state.mass^2

"""Radial problem in units of the current black-hole mass for a given background state."""
function current_problem(cache::EvolutionCache, state::BackgroundState)
    with_background(cache.template, KerrBackground(dimensionless_spin(state)), cache.initial_mass_coupling * state.mass)
end

"""Charges of a current-unit solution converted to reference units."""
reference_charges(charges::Charges, state::BackgroundState) =
    Charges(state.mass * charges.energy, state.mass^2 * charges.angular_momentum)

"""Losses per unit `τ` in reference units from current-unit fluxes."""
reference_losses(losses::Charges, state::BackgroundState) =
    Charges(losses.energy, state.mass * losses.angular_momentum)

"""Wave actions in reference units: `N_ref = b² N` (like `J`)."""
reference_actions(actions::AbstractVector, state::BackgroundState) = state.mass^2 .* actions

"""Real mode frequencies in reference units, `ω_ref = ω / b`."""
reference_frequencies(frequencies::AbstractVector, state::BackgroundState) = real.(frequencies) ./ state.mass

"""Real mode frequencies decoded from the unknown vector (current units)."""
real_mode_frequencies(layout::UnknownLayout, unknowns::AbstractVector) =
    real.(_complexify(unknowns[parameter_range(layout)][frequency_range(layout)]))

"""
Row ranges `(charges, actions, frequencies)` of `reference_observables`; `charges`
is empty without `include_charges`.
"""
function observable_ranges(mode_count::Int, include_charges::Bool)
    first_action = include_charges ? 3 : 1
    return (charges=1:(include_charges ? 2 : 0), actions=first_action:first_action + mode_count - 1,
            frequencies=first_action + mode_count:first_action + 2mode_count - 1)
end

"""
Everything the balances differentiate, stacked as one vector in reference
units: `[E, J, N₁ … N_K, ω₁ … ω_K]`. Without `include_charges` the energy and
angular momentum are not evaluated and the vector is `[N₁ … N_K, ω₁ … ω_K]`.
One dual pass of this function gives all derivatives at once.
"""
function reference_observables(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale,
                               state::BackgroundState; include_charges::Bool=true)
    totals = integrated_densities(problem, layout, unknowns, scale; charges=include_charges)
    component_rows = include_charges ? (3:length(totals)) : (1:length(totals))
    actions = reference_actions(phase_weights(problem.mode_set) * totals[component_rows], state)
    frequencies = reference_frequencies(real_mode_frequencies(layout, unknowns), state)
    include_charges || return vcat(actions, frequencies)
    charges = reference_charges(Charges(totals[1], totals[2]), state)
    return vcat(charges.energy, charges.angular_momentum, actions, frequencies)
end

"""Outer horizon radius `M + √(M² − a²)` of a hole with the given mass and angular momentum."""
outer_horizon_radius(mass, angular_momentum) = mass + sqrt(mass^2 - (angular_momentum / mass)^2)

"""
    inverse_metric_components(mass, angular_momentum, horizon_offset, x)

Boyer–Lindquist inverse Kerr metric for a hole of the given mass and angular
momentum (any consistent units), at radius `r = r₊ + horizon_offset` and
`x = cos θ`, as the vector `[g^{tt}, g^{tφ}, g^{φφ}, g^{rr}, g^{θθ}, ln ρ²]`.
The last entry gives `ln(−g) = 2 ln ρ² + ln sin²θ`. Taking the horizon offset
rather than `r` keeps `Δ = s(s + r₊ − r₋)` exact down the throat.
"""
function inverse_metric_components(mass, angular_momentum, horizon_offset, x)
    spin_length = angular_momentum / mass
    root = sqrt(mass^2 - spin_length^2)
    r = mass + root + horizon_offset
    sin2 = 1 - x^2
    horizon_function = horizon_offset * (horizon_offset + 2root)
    rho_squared = r^2 + spin_length^2 * x^2
    sum_r2a2 = r^2 + spin_length^2
    return [-(sum_r2a2^2 - spin_length^2 * horizon_function * sin2) / (rho_squared * horizon_function),
            -2mass * spin_length * r / (rho_squared * horizon_function),
            (horizon_function - spin_length^2 * sin2) / (rho_squared * horizon_function * sin2),
            horizon_function / rho_squared,
            1 / rho_squared,
            log(rho_squared)]
end

"""
    boyer_lindquist_metric_work(solution, state, state_rate)

Superseded by the action-based `metric_work`; kept for comparison. Rate of
work done on the cloud by the time-dependent metric, in reference units per
unit `τ`,

    w = ½ ∫ [ (∂_τ g^{αβ}) ∂_αΦ ∂_βΦ + (½ g^{αβ}∂_αΦ∂_βΦ + V) ∂_τ ln(−g) ] √−g d³x,

with `∂_τ` taken at fixed Boyer–Lindquist coordinates in reference units along
the direction `(b', j')` given by `state_rate`. This chart is singular at a
moving horizon: the integrand grows like `1/(r − r₊)`, so the result depends on
the inner endpoint (`validation/metric_work_endpoint.jl`).
"""
function boyer_lindquist_metric_work(solution::RadialSolution, state::BackgroundState, state_rate::BackgroundState)
    problem, layout, scale = solution.problem, solution.layout, solution.scale
    mesh, grid = problem.mesh, problem.quadrature.angular
    mass = state.mass
    initial_mass_coupling = problem.mass_coupling / mass
    frequencies = real_component_frequencies(problem, layout, solution.unknowns)
    weights = simpson_weights(mesh.rstar)

    total = 0.0
    for node in 1:layout.n_nodes
        offset = mesh.offsets[node]
        current_state = physical_state(layout, solution.unknowns, scale, node)
        bilinears = field_bilinears(problem, offset, current_state.profiles, current_state.slopes, frequencies)
        physical_offset = mass * offset
        shell = 0.0
        for (polar_index, x) in enumerate(grid.x)
            metric = inverse_metric_components(mass, state.angular_momentum, physical_offset, x)
            # ∂_τ at fixed Boyer–Lindquist r: the offset from the moving horizon changes by −∂_τ r₊,
            # formed so that the value part stays exactly `physical_offset`.
            metric_rate = ForwardDiff.derivative(0.0) do epsilon
                perturbed_mass = mass + epsilon * state_rate.mass
                perturbed_angular_momentum = state.angular_momentum + epsilon * state_rate.angular_momentum
                perturbed_offset = physical_offset + (outer_horizon_radius(mass, state.angular_momentum) -
                                                      outer_horizon_radius(perturbed_mass, perturbed_angular_momentum))
                inverse_metric_components(perturbed_mass, perturbed_angular_momentum, perturbed_offset, x)
            end
            # ⟨∂_αΦ ∂_βΦ⟩ in reference units: ∂_t = ∂_t̂ / b, ∂_r = ∂_r̂ / b.
            bilinear_vector = [bilinears.dt_dt[polar_index] / mass^2,
                               2bilinears.dt_dphi[polar_index] / mass,
                               bilinears.dphi_dphi[polar_index],
                               bilinears.dr_dr[polar_index] / mass^2,
                               bilinears.dtheta_dtheta[polar_index]]
            kinetic = dot(metric[1:5], bilinear_vector)
            lagrangian = kinetic / 2 + initial_mass_coupling^2 * bilinears.potential[polar_index]
            rho_squared = exp(metric[6])
            shell += grid.weights[polar_index] * rho_squared *
                     (dot(metric_rate[1:5], bilinear_vector) + lagrangian * 2metric_rate[6]) / 2
        end
        # √−g d³x = ρ² dr dx dφ with dr = b dr̂ = b (dr̂/dr*) dr*.
        total += weights[node] * mass * dr_drstar(problem.background, offset) * 2pi * shell
    end
    return total
end

"""
    BackgroundDerivatives

Total derivatives, at fixed amplitudes and fixed tortoise mesh in units of the
current mass, of the reference-unit charges, wave actions and real mode
frequencies with respect to the black-hole mass `b` and angular momentum `j`.
"""
struct BackgroundDerivatives{T<:Real}
    charges_wrt_mass::Charges{T}
    charges_wrt_angular_momentum::Charges{T}
    actions_wrt_mass::Vector{T}
    actions_wrt_angular_momentum::Vector{T}
    frequencies_wrt_mass::Vector{T}
    frequencies_wrt_angular_momentum::Vector{T}
end

"""
Split a stacked derivative of `reference_observables` along `b` and along `j`.
Without `include_charges` the charge derivatives are `NaN`.
"""
function BackgroundDerivatives(wrt_mass::AbstractVector, wrt_angular_momentum::AbstractVector, mode_count::Int;
                               include_charges::Bool=true)
    ranges = observable_ranges(mode_count, include_charges)
    charges_of(derivative) = include_charges ? Charges(derivative[1], derivative[2]) : Charges(NaN, NaN)
    BackgroundDerivatives(charges_of(wrt_mass), charges_of(wrt_angular_momentum),
                          wrt_mass[ranges.actions], wrt_angular_momentum[ranges.actions],
                          wrt_mass[ranges.frequencies], wrt_angular_momentum[ranges.frequencies])
end

"""
    background_tangents(cache, solution, state)

Matrix `[du/db du/dj]` (one column per background parameter) of the response
of the discretised solution to the black-hole mass and angular momentum at
fixed amplitudes and fixed tortoise mesh in units of the current mass, by the
implicit-function theorem `du/db = −J⁻¹ ∂F/∂b|_u` with the factorised Jacobian
of the converged solution. The partial `∂F/∂b|_u` is exact (ForwardDiff
through the background, mass coupling and mesh offsets: the outer-boundary
rows vary as `e^{−κr}` with `d(κr)/dα ≈ 10³`, far too steeply for finite
differences).
"""
function background_tangents(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState)
    layout, unknowns, scale, amplitudes = solution.layout, solution.unknowns, solution.scale, solution.amplitudes
    residual_partials = ForwardDiff.jacobian([0.0, 0.0]) do epsilon
        perturbed = BackgroundState(state.mass + epsilon[1], state.angular_momentum + epsilon[2])
        residual(current_problem(cache, perturbed), layout, unknowns, amplitudes, scale)
    end
    return -(solution.jacobian_factorization \ residual_partials)
end

"""Dual-number tag of `observable_partials` (all partials of one pass share it)."""
struct ObservableTag end

"""
    observable_partials(cache, solution, state, tangents, amplitude_directions; include_charges=true)

Value vector and derivative matrix of `reference_observables` from one dual
pass with `2 + K` partial slots: slot 1 is the black-hole mass `b`, slot 2 the
angular momentum `j`, and slot `2 + i` the amplitude `A_i`. The problem and
the state carry the explicit `b`, `j` dependence (seeded in slots 1 and 2) and
the unknown vector is displaced along the columns of `[tangents
amplitude_directions]`, so that the derivative along `b` and `j` is the total
derivative `∂O/∂b|_u + (∂O/∂u) du/db` and the derivative along `A_i` is
`(∂O/∂u) du/dA_i`. Returns `(values, partials)` with `partials[row, slot]`.
"""
function observable_partials(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState,
                             tangents::AbstractMatrix, amplitude_directions::AbstractMatrix; include_charges::Bool=true)
    directions = hcat(tangents, amplitude_directions)
    return _observable_partials(cache, solution, state, directions, Val(size(directions, 2)); include_charges)
end

function _observable_partials(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState,
                              directions::AbstractMatrix, ::Val{N}; include_charges::Bool) where {N}
    seed(value, slot) = ForwardDiff.Dual{ObservableTag}(value, ForwardDiff.Partials(ntuple(k -> k == slot ? 1.0 : 0.0, Val(N))))
    dual_state = BackgroundState(seed(state.mass, 1), seed(state.angular_momentum, 2))
    dual_unknowns = [ForwardDiff.Dual{ObservableTag}(value, ForwardDiff.Partials(ntuple(slot -> directions[index, slot], Val(N))))
                     for (index, value) in enumerate(solution.unknowns)]
    observables = reference_observables(current_problem(cache, dual_state), solution.layout, dual_unknowns, solution.scale,
                                        dual_state; include_charges)
    return ForwardDiff.value.(observables), [ForwardDiff.partials(observable, slot) for observable in observables, slot in 1:N]
end

"""
`SolutionTangents` of a converged solution on the background `state`: `amplitude_response` and
`background_tangents` (one dual residual pass).
"""
function solution_tangents(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState)
    background_directions = background_tangents(cache, solution, state)
    return SolutionTangents(solution, state.mass, state.angular_momentum, amplitude_response(solution),
                            background_directions[:, 1], background_directions[:, 2])
end

"""
    predicted_change(tangents, amplitudes, state)

First-order change of `tangents.solution`'s unknowns (in its units) for moving to `amplitudes` and background
`state`; pass it as `guess_change` to `solve_radial` with `guess = tangents.solution`.
"""
predicted_change(tangents::SolutionTangents, amplitudes::AbstractVector, state::BackgroundState) =
    tangents.wrt_amplitudes * (amplitudes .- tangents.solution.amplitudes) .+
    tangents.wrt_mass .* (state.mass - tangents.black_hole_mass) .+
    tangents.wrt_angular_momentum .* (state.angular_momentum - tangents.black_hole_angular_momentum)

"""
    background_charge_derivatives(cache, solution, state; relative_step=nothing, include_charges=true)

`BackgroundDerivatives` at fixed amplitudes and fixed tortoise mesh (in units
of the current mass). By the implicit-function theorem, for every observable `O`
of `reference_observables`,

    dO/db = ∂O/∂b|_u + (∂O/∂u) du/db,    du/db = −J⁻¹ ∂F/∂b|_u,

where `du/db` is `background_tangents` and both terms come from one dual pass
of the observables (`observable_partials`), exactly. With a number for
`relative_step` the explicit partial `∂O/∂b|_u` is instead a central
difference over the background with that relative step (the earlier
implementation, kept as a reference; benign for the observables, unlike for
the residual) and only the last term is a dual pass. Without
`include_charges` the energy and angular momentum are skipped and their
derivatives are `NaN`.
"""
function background_charge_derivatives(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                                       relative_step=nothing, include_charges::Bool=true)
    relative_step === nothing || return central_difference_background_derivatives(cache, solution, state, relative_step;
                                                                                  include_charges)
    tangents = background_tangents(cache, solution, state)
    _, partials = observable_partials(cache, solution, state, tangents, zeros(n_unknowns(solution.layout), 0); include_charges)
    return BackgroundDerivatives(partials[:, 1], partials[:, 2], n_modes(solution.problem); include_charges)
end

"""
    central_difference_background_derivatives(cache, solution, state, relative_step; include_charges=true)

Reference for `background_charge_derivatives`: the explicit partial of the
observables over the background is a central difference of
`reference_observables` at `±step` in `b` and in `j`; `via unknowns` is a dual
pass along each tangent.
"""
function central_difference_background_derivatives(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState,
                                                   relative_step; include_charges::Bool=true)
    problem, layout, unknowns, scale = solution.problem, solution.layout, solution.unknowns, solution.scale
    tangents = background_tangents(cache, solution, state)
    function derivative_along(slot::Int, direction::BackgroundState, step)
        plus = BackgroundState(state.mass + step * direction.mass, state.angular_momentum + step * direction.angular_momentum)
        minus = BackgroundState(state.mass - step * direction.mass, state.angular_momentum - step * direction.angular_momentum)
        explicit_partial = (reference_observables(current_problem(cache, plus), layout, unknowns, scale, plus; include_charges) .-
                            reference_observables(current_problem(cache, minus), layout, unknowns, scale, minus; include_charges)) ./ (2step)
        via_unknowns = ForwardDiff.derivative(0.0) do epsilon
            reference_observables(problem, layout, unknowns .+ epsilon .* tangents[:, slot], scale, state; include_charges)
        end
        return explicit_partial .+ via_unknowns
    end
    mass_step = relative_step * state.mass
    spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
    BackgroundDerivatives(derivative_along(1, BackgroundState(1.0, 0.0), mass_step),
                          derivative_along(2, BackgroundState(0.0, 1.0), spin_step), n_modes(problem); include_charges)
end

"""
    metric_work(derivatives, mode_frequencies, state_rate)

Rate of work done on the cloud by the changing black hole, in reference units
per unit `τ`, from the adiabatic invariance of the wave actions:

    w = (∂_b E − Σ_i ω_i ∂_b N_i) b' + (∂_j E − Σ_i ω_i ∂_j N_i) j' = (∂_b E)|_N b' + (∂_j E)|_N j',

with all derivatives at fixed amplitudes (`background_charge_derivatives`) and
`ω_i` the real reference-unit mode frequencies. The angular-momentum
counterpart vanishes identically because `J = Σ m_i N_i`. In the linear regime
(`E_i = ω_i N_i`) this reduces to `w = Σ_i N_i dω_i/dτ`. Unlike the
Boyer–Lindquist expression (`boyer_lindquist_metric_work`), no metric
derivative near the moving horizon is needed: the near-horizon densities of `E`
and `Σ ω_i N_i` coincide.
"""
function metric_work(derivatives::BackgroundDerivatives, mode_frequencies::AbstractVector, state_rate::BackgroundState)
    work_wrt_mass = derivatives.charges_wrt_mass.energy - sum(mode_frequencies .* derivatives.actions_wrt_mass)
    work_wrt_spin = derivatives.charges_wrt_angular_momentum.energy - sum(mode_frequencies .* derivatives.actions_wrt_angular_momentum)
    return work_wrt_mass * state_rate.mass + work_wrt_spin * state_rate.angular_momentum
end

"""
    BalanceTerms

Everything entering the evolving-background balances at one instant, in
reference units: amplitude derivatives of the charges and of the wave actions,
background derivatives at fixed amplitudes, charge and action losses, the
action-based metric work, and the background rates driven by the horizon fluxes.
Also the cloud's charges and actions, `background_directions`, the matrix
`[du/db du/dj]` of `background_tangents`, and `amplitude_directions`, the
matrix `[du/dA_i]` of `amplitude_response`. When the energy, angular
momentum and metric work were not requested (`_computes_charges`),
`charge_jacobian`, the charge derivatives and `metric_work` are `NaN`.
"""
struct BalanceTerms{T<:Real}
    charge_jacobian::Matrix{T}
    action_jacobian::Matrix{T}
    derivatives::BackgroundDerivatives{T}
    losses::Charges{T}
    action_losses::Vector{T}
    metric_work::T
    state_rate::BackgroundState{T}
    charges::Charges{T}
    actions::Vector{T}
    background_directions::Matrix{T}
    amplitude_directions::Matrix{T}
end

"""
Derivatives by one dual pass over the observables (`observable_partials`): the
`(derivatives, charge_matrix, action_matrix, charges, actions)` of `balance_terms`.
"""
function _exact_balance_pieces(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState,
                               tangents::AbstractMatrix, directions::AbstractMatrix, include_charges::Bool)
    mode_count = n_modes(solution.problem)
    values, partials = observable_partials(cache, solution, state, tangents, directions; include_charges)
    rows = observable_ranges(mode_count, include_charges)
    amplitude_slots = 3:2 + mode_count
    charges = include_charges ? Charges(values[1], values[2]) : reference_charges(cloud_charges(solution), state)
    return (derivatives=BackgroundDerivatives(partials[:, 1], partials[:, 2], mode_count; include_charges),
            charge_matrix=include_charges ? partials[rows.charges, amplitude_slots] : fill(NaN, 2, mode_count),
            action_matrix=partials[rows.actions, amplitude_slots], charges, actions=values[rows.actions])
end

"""
The same pieces as `_exact_balance_pieces`, with the background derivatives from
central finite differences over four re-solved radial problems (the slow
reference used for validation) and the Jacobians from `charge_jacobian` and `action_jacobian`.
"""
function _resolved_balance_pieces(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState,
                                  directions::AbstractMatrix, relative_step, include_charges::Bool)
    mode_count = n_modes(solution.problem)
    function observables_at(perturbed::BackgroundState)
        perturbed_solution = solve_radial(current_problem(cache, perturbed), solution.amplitudes; guess=solution, tol=cache.newton_tol)
        reference_observables(perturbed_solution.problem, perturbed_solution.layout, perturbed_solution.unknowns,
                              perturbed_solution.scale, perturbed; include_charges)
    end
    mass_step = relative_step * state.mass
    spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
    wrt_mass = (observables_at(BackgroundState(state.mass + mass_step, state.angular_momentum)) .-
                observables_at(BackgroundState(state.mass - mass_step, state.angular_momentum))) ./ (2mass_step)
    wrt_spin = (observables_at(BackgroundState(state.mass, state.angular_momentum + spin_step)) .-
                observables_at(BackgroundState(state.mass, state.angular_momentum - spin_step))) ./ (2spin_step)
    return (derivatives=BackgroundDerivatives(wrt_mass, wrt_spin, mode_count; include_charges),
            charge_matrix=include_charges ? charge_jacobian(solution; directions) .* [state.mass, state.mass^2] : fill(NaN, 2, mode_count),
            action_matrix=state.mass^2 .* action_jacobian(solution; directions),
            charges=reference_charges(cloud_charges(solution), state),
            actions=reference_actions(phase_actions(solution), state))
end

"""
    balance_terms(cache, solution, state; resolve_background=false, relative_step=1e-4,
                  include_charges=_computes_charges(cache))

Assemble the balance terms. The background rates are the horizon fluxes times
`(f_a/M_Pl)²` (see the module docstring). Background derivatives at fixed
amplitudes and tortoise mesh in units of the current mass, and the amplitude
Jacobians, come from one dual pass over the observables with the tangent
directions `du/db`, `du/dj` and the amplitude directions `du/dA_i`
(`observable_partials`); with `resolve_background=true` the background
derivatives are instead central finite differences over four re-solved radial
problems (the slow reference used for validation). The energy, the angular
momentum and everything derived from them (`charge_jacobian`, the charge
derivatives, `metric_work`) are computed only with `include_charges`, which by
default follows the cache (the energy balance, `include_metric_work` or
`record_metric_work`); the action balance does not use them.
"""
function balance_terms(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                       resolve_background=false, relative_step=1e-4, include_charges::Bool=_computes_charges(cache))
    inner_flux = outward_flux(solution, 1)
    state_rate = BackgroundState(-cache.coupling_to_hole * inner_flux.energy,
                                 -cache.coupling_to_hole * state.mass * inner_flux.angular_momentum)
    directions = amplitude_response(solution)
    tangents = background_tangents(cache, solution, state)
    pieces = resolve_background ? _resolved_balance_pieces(cache, solution, state, directions, relative_step, include_charges) :
                                  _exact_balance_pieces(cache, solution, state, tangents, directions, include_charges)
    frequencies = reference_frequencies(mode_frequencies(solution), state)
    work = include_charges ? metric_work(pieces.derivatives, frequencies, state_rate) : NaN
    BalanceTerms(pieces.charge_matrix, pieces.action_matrix, pieces.derivatives,
                 reference_losses(boundary_losses(solution), state), state.mass .* action_losses(solution), work, state_rate,
                 pieces.charges, pieces.actions, tangents, directions)
end

"""Tangents of `solution` from the directions its balance terms already hold (no further solves)."""
SolutionTangents(solution::RadialSolution, state::BackgroundState, terms::BalanceTerms) =
    SolutionTangents(solution, state.mass, state.angular_momentum, terms.amplitude_directions,
                     terms.background_directions[:, 1], terms.background_directions[:, 2])

"""`SolutionTangents` from the balance terms when the cache uses the tangent predictor, otherwise `nothing`."""
_tangents_from_terms(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState, terms::BalanceTerms) =
    cache.tangent_predictor ? SolutionTangents(solution, state, terms) : nothing

"""
    amplitude_rates(terms::BalanceTerms; balance=:action, include_metric_work=false)

`A'` from the energy balance
`[∂_A E; ∂_A J] A' = −L + (w, 0) − (∂_b Q b' + ∂_j Q j')` (`w` only with
`include_metric_work`), or from the action balance
`[∂_A N] A' = −L_N − (∂_b N b' + ∂_j N j')`, which needs no work term.
"""
function amplitude_rates(terms::BalanceTerms; balance::Symbol=:action, include_metric_work::Bool=false)
    derivatives, rate = terms.derivatives, terms.state_rate
    if balance == :action
        drift = rate.mass .* derivatives.actions_wrt_mass .+ rate.angular_momentum .* derivatives.actions_wrt_angular_momentum
        return terms.action_jacobian \ (-terms.action_losses .- drift)
    end
    balance == :energy || throw(ArgumentError("balance must be :energy or :action; got :$balance"))
    any(isnan, terms.charge_jacobian) &&
        throw(ArgumentError("these balance terms were assembled without the energy and angular momentum; build the cache " *
                            "with record_metric_work=true (or balance=:energy) to use the energy balance"))
    drift = rate.mass * derivatives.charges_wrt_mass + rate.angular_momentum * derivatives.charges_wrt_angular_momentum
    work = include_metric_work ? terms.metric_work : 0.0
    return terms.charge_jacobian \ (-collect(terms.losses) .+ [work, 0.0] .- collect(drift))
end

"""`amplitude_rates` with the balance and metric-work choice of an evolution cache."""
amplitude_rates(terms::BalanceTerms, cache::EvolutionCache) =
    amplitude_rates(terms; cache.balance, include_metric_work=cache.include_metric_work)

"""
Right-hand side of the spindown evolution. Every radial solve is warm-started
from the solution at the last accepted step (moved by `predicted_change` when the
cache has its tangents), never from a stage of a rejected step: after a rejection the retry's stages can be far from the rejected trial
stages in both amplitudes and background, and that combination broke the warm
start (`validation/spindown_warmstart_check.jl`). If a solve still fails, the
rates are set to zero, the remaining stages of the step are skipped, and
`_reject_failed_step` makes the integrator retry with a smaller step.
"""
function _spindown_rhs!(rates, state_vector, cache::EvolutionCache, time)
    if cache.solve_failed
        rates .= 0
        return nothing
    end
    evaluation = _spindown_evaluation(state_vector, cache, cache.accepted_solution, cache.accepted_tangents)
    cache.solve_counts.stage += 1
    solution, state, terms = evaluation.solution, evaluation.state, evaluation.terms
    if evaluation.rates === nothing
        cache.solve_counts.failed += 1
        cache.solve_failed = true
        rates .= 0
        cache.verbose && (println("t = ", time, "  radial solve failed (|F| = ", solution.residual_norm,
                                  "); the step will be rejected"); flush(stdout))
        return nothing
    end
    rates .= evaluation.rates
    cache.last_solution, cache.last_rates, cache.last_tangents = solution, evaluation.rates, evaluation.tangents
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(evaluation.amplitudes), mode_frequencies(solution),
                                                  terms.charges, terms.losses, state.mass, dimensionless_spin(state),
                                                  terms.metric_work, terms.actions, terms.action_losses))
    _report_step(cache, time, evaluation.amplitudes, solution,
                 "  b = $(state.mass)  j = $(state.angular_momentum)  w = $(terms.metric_work)")
    return nothing
end

"""Integrator-state rates `[(ln) A' …, b', j']` from the balance terms at `amplitudes`."""
_spindown_state_rates(terms::BalanceTerms, amplitudes, cache::EvolutionCache) =
    vcat(_state_rates(amplitude_rates(terms, cache), amplitudes, cache.log_amplitudes),
         terms.state_rate.mass, terms.state_rate.angular_momentum)

"""
Radial solve warm-started from `guess` and balance terms at the integrator
state `state_vector`. Returns the amplitudes, background state, radial
solution, balance terms, integrator-state rates and the solution's tangents;
`rates`, `terms` and `tangents` are `nothing` if the solve did not converge.
With the tangent predictor and `guess_tangents` belonging to `guess`, the solve
starts from `guess` moved by `predicted_change`. Mutates nothing, so evaluations
can run concurrently on one cache.
"""
function _spindown_evaluation(state_vector, cache::EvolutionCache, guess, guess_tangents=nothing)
    mode_count = n_modes(cache.template)
    amplitudes = _state_amplitudes(state_vector, mode_count, cache.log_amplitudes)
    state = BackgroundState(state_vector[mode_count + 1], state_vector[mode_count + 2])
    predict = cache.tangent_predictor && guess_tangents !== nothing && guess_tangents.solution === guess
    guess_change = predict ? predicted_change(guess_tangents, amplitudes, state) : nothing
    solution = solve_radial(current_problem(cache, state), amplitudes; guess, guess_change, tol=cache.newton_tol)
    solution.converged || return (; amplitudes, state, solution, terms=nothing, rates=nothing, tangents=nothing)
    terms = balance_terms(cache, solution, state)
    return (; amplitudes, state, solution, terms, rates=_spindown_state_rates(terms, amplitudes, cache),
            tangents=_tangents_from_terms(cache, solution, state, terms))
end

"""
Evaluation at `base_state + perturbation` for a Jacobian column, solved from
`base_solution`, the converged solution at `base_state`, moved by its tangents
(`cache.accepted_tangents`, set with the base by `_accepted_evaluation!`). All
perturbed solves go through this function.
"""
_perturbed_spindown_evaluation(cache::EvolutionCache, base_solution::RadialSolution, base_state, perturbation) =
    _spindown_evaluation(base_state .+ perturbation, cache, base_solution, cache.accepted_tangents)

"""
Radial solution and rates at the last accepted state `state_vector`. After a
step of a first-same-as-last integrator the cache already holds them; otherwise
it holds the last stage of the accepted step, from which one solve at the
accepted state is done (and kept). If that solve fails the evolution stops with
an error, since an accepted step cannot be rejected afterwards.
"""
function _accepted_evaluation!(cache::EvolutionCache, state_vector, time)
    _solution_at_state(cache.accepted_solution, state_vector, cache) &&
        return cache.accepted_solution, cache.accepted_rates
    evaluation = _spindown_evaluation(state_vector, cache, cache.accepted_solution, cache.accepted_tangents)
    cache.solve_counts.accepted_state += 1
    if evaluation.rates === nothing
        cache.solve_counts.failed += 1
        error("radial solve at the accepted state t = $time did not converge")
    end
    cache.accepted_solution, cache.accepted_rates, cache.accepted_tangents =
        evaluation.solution, evaluation.rates, evaluation.tangents
    _report_line(cache, "accepted t = $time: radial solve at the accepted state, newton = $(evaluation.solution.iterations)")
    return cache.accepted_solution, cache.accepted_rates
end

"""
Jacobian request of the spindown integrator, which integrators make at the last
accepted state; its solution and rates (`_accepted_evaluation!`) are the base.
"""
function _spindown_jacobian!(jacobian_matrix, state_vector, cache::EvolutionCache, settings::FiniteDifferenceJacobian,
                             store::JacobianStore, time)
    base() = _accepted_evaluation!(cache, state_vector, time)
    _fill_jacobian!(jacobian_matrix, state_vector, time, cache, settings, store, base,
                    (base_solution, perturbation) ->
                        _perturbed_spindown_evaluation(cache, base_solution, state_vector, perturbation))
end

"""
`isoutofdomain` for the spindown integrator: reject the step if any radial
solve in it failed (the integrator then shrinks the step by `qmin`), and clear
the flag for the next attempt. Called once per step attempt, after the last stage.
"""
function _reject_failed_step(state_vector, cache::EvolutionCache, time)
    failed = cache.solve_failed
    cache.solve_failed = false
    return failed
end

"""
Callback after every accepted step: keep the step's last stage solution, its
rates and its tangents as the warm start of the next step's solves. Stages of
rejected attempts are never used, since after a rejection the retry's stages can
be far from the rejected ones (`_spindown_rhs!`). A first-same-as-last integrator
(Tsit5, Rosenbrock23) evaluates the right-hand side at the accepted state as its
last stage, so this is the solution there; for other integrators (Rodas5P, FBDF)
it is a nearby stage, and the solve at the accepted state is done only when a
Jacobian needs it (`_accepted_evaluation!`).
"""
function _record_accepted_solution!(integrator)
    cache = integrator.p
    cache.accepted_solution, cache.accepted_rates, cache.accepted_tangents =
        cache.last_solution, cache.last_rates, cache.last_tangents
    u_modified!(integrator, false)
    return nothing
end

"""
    evolve_with_spindown(problem, amplitudes, time_span, decay_constant_ratio; balance=:action,
                         include_metric_work=false, record_metric_work=false, log_amplitudes=true, initial_step_fraction=0.05,
                         newton_tol=1e-10, verbose=false, initial_background, initial_guess=nothing,
                         tangent_predictor=true, algorithm=Tsit5(), jacobian_settings=FiniteDifferenceJacobian(),
                         callback=nothing, reltol, abstol, kwargs...)

Integrate `[A₁ … A_K, b, j]` (or `[ln A₁ … ln A_K, b, j]`, the default) with
`algorithm` (default Tsit5), starting from `initial_background`, by default `b = 1` and `j = a` of
`problem.background`. `problem`
defines the reference units: its `mass_coupling` is `α₀ = M₀ μ` and its mesh is
reused in units of the current mass at every step. `decay_constant_ratio =
f_a / M_Pl` (non-reduced Planck mass) sets how strongly the cloud's horizon
fluxes change the hole. Tolerances, the initial step and `log_amplitudes` are
as in `evolve_fixed_background`.

The default `balance=:action` uses the wave-action balance (any number of
modes, no metric-work term). With `balance=:energy`, `include_metric_work` adds
the action-based `metric_work` to the energy equation; it is off by default so
that earlier energy-balance runs reproduce, and is of relative order `(f_a/M_Pl)² E` like the
background-drift terms. The energy, the angular momentum and the work are computed (and the work
recorded in the diagnostics) only for the energy balance, with `include_metric_work`, or with
`record_metric_work=true`; under the default action balance the diagnostics hold `NaN` for the work
and the amplitude rates are unchanged (`_computes_charges`).

Each radial solve is warm-started from the solution at the last accepted step;
with `tangent_predictor` (default) it starts from that solution moved to first
order in the changes of amplitudes, `b` and `j` (`SolutionTangents`,
`predicted_change`), which costs one `solution_tangents` per accepted step and
replaces the Newton iterations a plain warm start needs for the change of
background. A step in which a solve fails is rejected and retried with a smaller step
(`_spindown_rhs!`). If the step size then reaches the integrator's minimum, the
integration stops with a non-success return code instead of using an
unconverged solution. `isoutofdomain` is used internally. A `callback` (for
example one that writes checkpoints) runs after the accepted-step bookkeeping,
so `integrator.p.accepted_solution` is the radial solution at `integrator.u`.

On the |211⟩ plateau the hole relaxes onto the superradiance threshold on a
time scale far shorter than the slow evolution, and Tsit5's step is held at its
stability limit (`docs/session-2026-11-12.md`). A stiff integrator such as
`algorithm=Rosenbrock23()` (re-exported from `OrdinaryDiffEqRosenbrock`) then
takes steps set by accuracy. It receives the finite-difference Jacobian of
`jacobian_settings` (`FiniteDifferenceJacobian`) and a zero time gradient. The
diagnostics log holds the stage evaluations only; `result.solve_counts`
separates them from the Jacobian's perturbed solves and from the solves at
accepted states that integrators which are not first-same-as-last (Rodas5P)
need before a new Jacobian, and the verbose log prints the latter two on lines
starting with `jacobian` and `accepted`.

To restart from a saved state `(t, A, b, j)`, pass the original run's `problem`
(it fixes the reference units and the mesh), the saved amplitudes, `time_span =
(t, t_final)` and `initial_background = BackgroundState(b, j)`. `initial_guess`
is a converged `RadialSolution` at the starting state to begin the first solve
from (for example from amplitude continuation); without it the first solve
starts cold, which can fail in the nonlinear regime.
"""
function evolve_with_spindown(problem::RadialProblem, amplitudes::AbstractVector, time_span, decay_constant_ratio::Real;
                              balance::Symbol=:action, include_metric_work=false, record_metric_work=false, log_amplitudes=true,
                              initial_step_fraction=0.05, newton_tol=1e-10, verbose=false,
                              initial_background::BackgroundState=BackgroundState(1.0, problem.background.a),
                              initial_guess=nothing, tangent_predictor=true, algorithm=Tsit5(),
                              jacobian_settings::FiniteDifferenceJacobian=FiniteDifferenceJacobian(), callback=nothing,
                              reltol=(log_amplitudes ? 1e-7 : 1e-6), abstol=(log_amplitudes ? 1e-7 : 1e-12), kwargs...)
    cache = EvolutionCache(problem; decay_constant_ratio, include_metric_work, record_metric_work, newton_tol, verbose, log_amplitudes,
                           balance, last_solution=initial_guess, tangent_predictor)
    initial_amplitudes = collect(float.(amplitudes))
    initial_solution = _solve_cached!(cache, current_problem(cache, initial_background), initial_amplitudes)
    cache.solve_counts.accepted_state += 1
    initial_solution.converged || error("initial radial solve did not converge")
    initial_terms = balance_terms(cache, initial_solution, initial_background)
    initial_rates = amplitude_rates(initial_terms, cache)
    cache.accepted_solution = initial_solution
    cache.accepted_rates = cache.last_rates = _spindown_state_rates(initial_terms, initial_amplitudes, cache)
    cache.accepted_tangents = cache.last_tangents = _tangents_from_terms(cache, initial_solution, initial_background, initial_terms)
    solver_options = _with_initial_step(kwargs, _initial_step(initial_amplitudes, initial_rates, initial_step_fraction))
    initial_state = vcat(log_amplitudes ? log.(initial_amplitudes) : initial_amplitudes,
                         [initial_background.mass, initial_background.angular_momentum])
    store = JacobianStore(length(initial_state))
    jacobian! = (matrix, state, cache, time) -> _spindown_jacobian!(matrix, state, cache, jacobian_settings, store, time)
    ode = ODEProblem(_ode_function(_spindown_rhs!, jacobian!), initial_state, time_span, cache)
    accepted_step_callback = DiscreteCallback((state, time, integrator) -> true, _record_accepted_solution!;
                                              save_positions=(false, false))
    callbacks = callback === nothing ? accepted_step_callback : CallbackSet(accepted_step_callback, callback)
    ode_solution = solve(ode, algorithm; reltol, abstol, callback=callbacks, isoutofdomain=_reject_failed_step,
                         solver_options...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution, log_amplitudes, cache.solve_counts)
end
