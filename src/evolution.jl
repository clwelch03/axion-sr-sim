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
    EvolutionResult

Output of a driver: the `OrdinaryDiffEq` solution object (state `[A₁, …, A_K]` or
`[A₁, …, A_K, b, j]`, with `ln A` in place of `A` when `log_amplitudes` is true),
the diagnostics log, and the final radial solution. Read the amplitudes with
`amplitude_trajectory` rather than from `ode_solution.u` directly.
"""
struct EvolutionResult{S, T<:Real, R<:RadialSolution}
    ode_solution::S
    diagnostics::Vector{EvolutionDiagnostics{T}}
    final_radial_solution::R
    log_amplitudes::Bool
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

"""Mutable scratch shared by the right-hand-side evaluations of one evolution."""
mutable struct EvolutionCache{P<:RadialProblem, T<:Real}
    template::P
    initial_mass_coupling::T
    coupling_to_hole::T   # (f_a / M_Pl)²: converts cloud fluxes in f_a² units into black-hole mass rates
    include_metric_work::Bool
    last_solution::Union{Nothing, RadialSolution}
    diagnostics::Vector{EvolutionDiagnostics{T}}
    newton_tol::T
    verbose::Bool
    log_amplitudes::Bool
    balance::Symbol       # :energy or :action
    accepted_solution::Union{Nothing, RadialSolution}   # radial solution at the last accepted step (spindown warm start)
    solve_failed::Bool    # a radial solve in the current step attempt did not converge (spindown rejects the step)
end

"""
    EvolutionCache(problem; decay_constant_ratio=0.0, include_metric_work=false, newton_tol=1e-10,
                   verbose=false, log_amplitudes=true, balance=:action, last_solution=nothing)

Fresh cache for an evolution of `problem` (also handy for calling
`balance_terms` directly).
"""
function EvolutionCache(problem::RadialProblem; decay_constant_ratio=0.0, include_metric_work=false, newton_tol=1e-10,
                        verbose=false, log_amplitudes=true, balance=:action, last_solution=nothing)
    _check_balance(problem, balance)
    EvolutionCache(problem, problem.mass_coupling, float(decay_constant_ratio)^2, include_metric_work, last_solution,
                   EvolutionDiagnostics{Float64}[], float(newton_tol), verbose, log_amplitudes, balance,
                   last_solution, false)
end

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

"""Radial solve warm-started from the cache's previous solution."""
function _solve_cached!(cache::EvolutionCache, problem::RadialProblem, amplitudes)
    solution = solve_radial(problem, amplitudes; guess=cache.last_solution, tol=cache.newton_tol)
    cache.last_solution = solution
    return solution
end

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

function _fixed_background_rhs!(rates, state, cache::EvolutionCache, time)
    amplitudes = _state_amplitudes(state, n_modes(cache.template), cache.log_amplitudes)
    solution = _solve_cached!(cache, cache.template, amplitudes)
    rates .= _state_rates(amplitude_rates(solution; cache.balance), amplitudes, cache.log_amplitudes)
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(amplitudes), mode_frequencies(solution),
                                                  cloud_charges(solution), boundary_losses(solution),
                                                  1.0, cache.template.background.a, 0.0,
                                                  phase_actions(solution), action_losses(solution)))
    _report_step(cache, time, amplitudes, solution, "")
    return nothing
end

"""
    evolve_fixed_background(problem, amplitudes, time_span; balance=:action, log_amplitudes=true,
                            initial_step_fraction=0.05, newton_tol=1e-10, verbose=false, reltol, abstol, kwargs...)

Integrate the amplitudes with Tsit5 on a frozen Kerr background, with the
action (default, any number of modes) or energy (two modes) balance.
With `log_amplitudes` (default) the state is `ln A` and the default tolerances
`reltol = abstol = 1e-7` bound the per-step relative error of the amplitudes
by about `1e-7 (1 + |ln A|)`; otherwise the state is `A` with the old defaults
`reltol = 1e-6`, `abstol = 1e-12`. The first step is `initial_step_fraction`
times the shortest initial e-folding time unless `dt` is passed
(`initial_step_fraction=nothing` leaves it to Tsit5's heuristic).
`verbose` prints one line per right-hand-side evaluation. Other keyword
arguments are passed to `solve` (e.g. `saveat`, `dtmax`).
"""
function evolve_fixed_background(problem::RadialProblem, amplitudes::AbstractVector, time_span; balance::Symbol=:action,
                                 log_amplitudes=true, initial_step_fraction=0.05, newton_tol=1e-10, verbose=false,
                                 reltol=(log_amplitudes ? 1e-7 : 1e-6), abstol=(log_amplitudes ? 1e-7 : 1e-12), kwargs...)
    cache = EvolutionCache(problem; newton_tol, verbose, log_amplitudes, balance)
    initial_amplitudes = collect(float.(amplitudes))
    initial_rates = amplitude_rates(_solve_cached!(cache, problem, initial_amplitudes); balance)
    solver_options = _with_initial_step(kwargs, _initial_step(initial_amplitudes, initial_rates, initial_step_fraction))
    initial_state = log_amplitudes ? log.(initial_amplitudes) : initial_amplitudes
    ode = ODEProblem(_fixed_background_rhs!, initial_state, time_span, cache)
    ode_solution = solve(ode, Tsit5(); reltol, abstol, solver_options...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution, log_amplitudes)
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
Everything the balances differentiate, stacked as one vector in reference
units: `[E, J, N₁ … N_K, ω₁ … ω_K]`. One dual pass of this function gives all
background derivatives at once.
"""
function reference_observables(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale,
                               state::BackgroundState)
    charges = reference_charges(cloud_charges(problem, layout, unknowns, scale), state)
    vcat(charges.energy, charges.angular_momentum,
         reference_actions(phase_actions(problem, layout, unknowns, scale), state),
         reference_frequencies(real_mode_frequencies(layout, unknowns), state))
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

"""Split a stacked derivative of `reference_observables` along `b` and along `j`."""
function BackgroundDerivatives(wrt_mass::AbstractVector, wrt_angular_momentum::AbstractVector, mode_count::Int)
    actions = 3:2 + mode_count
    frequencies = 3 + mode_count:2 + 2mode_count
    BackgroundDerivatives(Charges(wrt_mass[1], wrt_mass[2]), Charges(wrt_angular_momentum[1], wrt_angular_momentum[2]),
                          wrt_mass[actions], wrt_angular_momentum[actions],
                          wrt_mass[frequencies], wrt_angular_momentum[frequencies])
end

"""
    background_charge_derivatives(cache, solution, state; relative_step=1e-6)

`BackgroundDerivatives` at fixed amplitudes and fixed tortoise mesh (in units
of the current mass). By the implicit-function theorem, for every observable `O`
of `reference_observables`,

    dO/db = ∂O/∂b|_u + (∂O/∂u) du/db,    du/db = −J⁻¹ ∂F/∂b|_u,

where `∂F/∂b|_u` is exact (ForwardDiff through the background, mass coupling
and mesh offsets: the outer-boundary rows vary as `e^{−κr}` with
`d(κr)/dα ≈ 10³`, far too steeply for finite differences), `∂O/∂b|_u` is a
central difference over the background (benign), `J` is the factorised
Jacobian of the converged solution, and the last term is one dual pass along
`du/db`, exactly as in `charge_jacobian`.
"""
function background_charge_derivatives(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                                       relative_step=1e-6)
    problem, layout, unknowns, scale, amplitudes = solution.problem, solution.layout, solution.unknowns, solution.scale,
                                                   solution.amplitudes
    function derivative_along(direction::BackgroundState, step)
        plus = BackgroundState(state.mass + step * direction.mass, state.angular_momentum + step * direction.angular_momentum)
        minus = BackgroundState(state.mass - step * direction.mass, state.angular_momentum - step * direction.angular_momentum)
        residual_partial = ForwardDiff.derivative(0.0) do epsilon
            perturbed = BackgroundState(state.mass + epsilon * direction.mass,
                                        state.angular_momentum + epsilon * direction.angular_momentum)
            residual(current_problem(cache, perturbed), layout, unknowns, amplitudes, scale)
        end
        explicit_partial = (reference_observables(current_problem(cache, plus), layout, unknowns, scale, plus) .-
                            reference_observables(current_problem(cache, minus), layout, unknowns, scale, minus)) ./ (2step)
        unknown_direction = -(solution.jacobian_factorization \ residual_partial)
        via_unknowns = ForwardDiff.derivative(0.0) do epsilon
            reference_observables(problem, layout, unknowns .+ epsilon .* unknown_direction, scale, state)
        end
        return explicit_partial .+ via_unknowns
    end
    mass_step = relative_step * state.mass
    spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
    BackgroundDerivatives(derivative_along(BackgroundState(1.0, 0.0), mass_step),
                          derivative_along(BackgroundState(0.0, 1.0), spin_step), n_modes(problem))
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
"""
struct BalanceTerms{T<:Real}
    charge_jacobian::Matrix{T}
    action_jacobian::Matrix{T}
    derivatives::BackgroundDerivatives{T}
    losses::Charges{T}
    action_losses::Vector{T}
    metric_work::T
    state_rate::BackgroundState{T}
end

"""
    balance_terms(cache, solution, state; resolve_background=false, relative_step=1e-4)

Assemble the balance terms. The background rates are the horizon fluxes times
`(f_a/M_Pl)²` (see the module docstring). Background derivatives at fixed
amplitudes and tortoise mesh in units of the current mass come from
`background_charge_derivatives`; with `resolve_background=true` they are instead
central finite differences over four re-solved radial problems (the slow
reference used for validation). The metric work (`metric_work`) is always
computed and reported, and enters the energy balance only when the cache says so.
"""
function balance_terms(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                       resolve_background=false, relative_step=1e-4)
    amplitudes, mode_count = solution.amplitudes, n_modes(solution.problem)
    inner_flux = outward_flux(solution, 1)
    state_rate = BackgroundState(-cache.coupling_to_hole * inner_flux.energy,
                                 -cache.coupling_to_hole * state.mass * inner_flux.angular_momentum)

    derivatives = if resolve_background
        function observables_at(perturbed::BackgroundState)
            perturbed_solution = solve_radial(current_problem(cache, perturbed), amplitudes; guess=solution, tol=cache.newton_tol)
            reference_observables(perturbed_solution.problem, perturbed_solution.layout, perturbed_solution.unknowns,
                                  perturbed_solution.scale, perturbed)
        end
        mass_step = relative_step * state.mass
        spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
        wrt_mass = (observables_at(BackgroundState(state.mass + mass_step, state.angular_momentum)) .-
                    observables_at(BackgroundState(state.mass - mass_step, state.angular_momentum))) ./ (2mass_step)
        wrt_spin = (observables_at(BackgroundState(state.mass, state.angular_momentum + spin_step)) .-
                    observables_at(BackgroundState(state.mass, state.angular_momentum - spin_step))) ./ (2spin_step)
        BackgroundDerivatives(wrt_mass, wrt_spin, mode_count)
    else
        background_charge_derivatives(cache, solution, state)
    end

    directions = amplitude_response(solution)
    charge_matrix = charge_jacobian(solution; directions) .* [state.mass, state.mass^2]
    action_matrix = state.mass^2 .* action_jacobian(solution; directions)
    frequencies = reference_frequencies(mode_frequencies(solution), state)
    BalanceTerms(charge_matrix, action_matrix, derivatives, reference_losses(boundary_losses(solution), state),
                 state.mass .* action_losses(solution), metric_work(derivatives, frequencies, state_rate), state_rate)
end

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
    drift = rate.mass * derivatives.charges_wrt_mass + rate.angular_momentum * derivatives.charges_wrt_angular_momentum
    work = include_metric_work ? terms.metric_work : 0.0
    return terms.charge_jacobian \ (-collect(terms.losses) .+ [work, 0.0] .- collect(drift))
end

"""`amplitude_rates` with the balance and metric-work choice of an evolution cache."""
amplitude_rates(terms::BalanceTerms, cache::EvolutionCache) =
    amplitude_rates(terms; cache.balance, include_metric_work=cache.include_metric_work)

"""
Right-hand side of the spindown evolution. Every radial solve is warm-started
from the solution at the last accepted step, never from a stage of a rejected
step: after a rejection the retry's stages can be far from the rejected trial
stages in both amplitudes and background, and that combination broke the warm
start (`validation/spindown_warmstart_check.jl`). If a solve still fails, the
rates are set to zero, the remaining stages of the step are skipped, and
`_reject_failed_step` makes the integrator retry with a smaller step.
"""
function _spindown_rhs!(rates, state_vector, cache::EvolutionCache, time)
    mode_count = n_modes(cache.template)
    if cache.solve_failed
        rates .= 0
        return nothing
    end
    amplitudes = _state_amplitudes(state_vector, mode_count, cache.log_amplitudes)
    state = BackgroundState(state_vector[mode_count + 1], state_vector[mode_count + 2])
    solution = solve_radial(current_problem(cache, state), amplitudes; guess=cache.accepted_solution, tol=cache.newton_tol)
    if !solution.converged
        cache.solve_failed = true
        rates .= 0
        cache.verbose && (println("t = ", time, "  radial solve failed (|F| = ", solution.residual_norm,
                                  "); the step will be rejected"); flush(stdout))
        return nothing
    end
    cache.last_solution = solution
    terms = balance_terms(cache, solution, state)
    rates[1:mode_count] = _state_rates(amplitude_rates(terms, cache), amplitudes, cache.log_amplitudes)
    rates[mode_count + 1] = terms.state_rate.mass
    rates[mode_count + 2] = terms.state_rate.angular_momentum
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(amplitudes), mode_frequencies(solution),
                                                  reference_charges(cloud_charges(solution), state), terms.losses,
                                                  state.mass, dimensionless_spin(state), terms.metric_work,
                                                  reference_actions(phase_actions(solution), state), terms.action_losses))
    _report_step(cache, time, amplitudes, solution, "  b = $(state.mass)  j = $(state.angular_momentum)  w = $(terms.metric_work)")
    return nothing
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
Callback after every accepted step: keep the step's radial solution as the
warm start for the next step. Tsit5 evaluates the right-hand side at the
accepted state as its last stage (first-same-as-last), so `cache.last_solution`
is the solution there; the check makes a different integrator fail loudly.
"""
function _record_accepted_solution!(integrator)
    cache = integrator.p
    amplitudes = _state_amplitudes(integrator.u, n_modes(cache.template), cache.log_amplitudes)
    isapprox(cache.last_solution.amplitudes, amplitudes; rtol=1e-12) ||
        error("last radial solution is not at the accepted state (amplitudes $(cache.last_solution.amplitudes) " *
              "vs $amplitudes); the accepted-step warm start assumes a first-same-as-last integrator such as Tsit5")
    cache.accepted_solution = cache.last_solution
    u_modified!(integrator, false)
    return nothing
end

"""
    evolve_with_spindown(problem, amplitudes, time_span, decay_constant_ratio; balance=:action,
                         include_metric_work=false, log_amplitudes=true, initial_step_fraction=0.05,
                         newton_tol=1e-10, verbose=false, initial_background, initial_guess=nothing,
                         reltol, abstol, kwargs...)

Integrate `[A₁ … A_K, b, j]` (or `[ln A₁ … ln A_K, b, j]`, the default) with
Tsit5, starting from `initial_background`, by default `b = 1` and `j = a` of
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
background-drift terms. The work is recorded in the diagnostics either way.

Each radial solve is warm-started from the solution at the last accepted step,
and a step in which a solve fails is rejected and retried with a smaller step
(`_spindown_rhs!`). If the step size then reaches the integrator's minimum, the
integration stops with a non-success return code instead of using an
unconverged solution. Do not pass `callback` or `isoutofdomain`; they are used.

To restart from a saved state `(t, A, b, j)`, pass the original run's `problem`
(it fixes the reference units and the mesh), the saved amplitudes, `time_span =
(t, t_final)` and `initial_background = BackgroundState(b, j)`. `initial_guess`
is a converged `RadialSolution` at the starting state to begin the first solve
from (for example from amplitude continuation); without it the first solve
starts cold, which can fail in the nonlinear regime.
"""
function evolve_with_spindown(problem::RadialProblem, amplitudes::AbstractVector, time_span, decay_constant_ratio::Real;
                              balance::Symbol=:action, include_metric_work=false, log_amplitudes=true,
                              initial_step_fraction=0.05, newton_tol=1e-10, verbose=false,
                              initial_background::BackgroundState=BackgroundState(1.0, problem.background.a),
                              initial_guess=nothing,
                              reltol=(log_amplitudes ? 1e-7 : 1e-6), abstol=(log_amplitudes ? 1e-7 : 1e-12), kwargs...)
    cache = EvolutionCache(problem; decay_constant_ratio, include_metric_work, newton_tol, verbose, log_amplitudes, balance,
                           last_solution=initial_guess)
    initial_amplitudes = collect(float.(amplitudes))
    initial_solution = _solve_cached!(cache, current_problem(cache, initial_background), initial_amplitudes)
    initial_solution.converged || error("initial radial solve did not converge")
    cache.accepted_solution = initial_solution
    initial_rates = amplitude_rates(balance_terms(cache, initial_solution, initial_background), cache)
    solver_options = _with_initial_step(kwargs, _initial_step(initial_amplitudes, initial_rates, initial_step_fraction))
    initial_state = vcat(log_amplitudes ? log.(initial_amplitudes) : initial_amplitudes,
                         [initial_background.mass, initial_background.angular_momentum])
    ode = ODEProblem(_spindown_rhs!, initial_state, time_span, cache)
    accepted_step_callback = DiscreteCallback((state, time, integrator) -> true, _record_accepted_solution!;
                                              save_positions=(false, false))
    ode_solution = solve(ode, Tsit5(); reltol, abstol, callback=accepted_step_callback,
                         isoutofdomain=_reject_failed_step, solver_options...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution, log_amplitudes)
end
