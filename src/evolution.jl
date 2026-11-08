"""
Adiabatic time evolution of the cloud amplitudes.

At every instant the cloud is the stationary multiphase solution at the current
amplitudes; energy and angular-momentum balance between the slowly changing
charges and the boundary fluxes then fixes the amplitude rates. Two levels are
provided:

* `evolve_fixed_background`: the black hole does not respond (toy model),
  `[∂_A E; ∂_A J] A' = −L`.
* `evolve_with_spindown`: the black-hole mass and spin absorb the horizon
  fluxes, and the balance includes the background derivatives of the charges
  and the rate of metric work on the cloud.

Reference units are the *initial* black-hole mass `M₀` (and `f_a`): with
`b = M/M₀`, `j = J_BH/M₀²`, `τ = t/M₀`, a radial solution in units of the
current mass carries `E_ref = b E`, `J_ref = b² J`, energy losses per unit `τ`
equal to the current-unit fluxes, and angular-momentum losses `b` times them.
Two independent conservation laws give two amplitude equations, so these
drivers require exactly two bound modes.

Cloud charges are measured in units of `f_a² M₀` while the black hole's are in
units of `M₀` (with `G = 1`), so the horizon fluxes change the hole's mass and
spin only after multiplication by `(f_a / M_Pl)²`, where `M_Pl = G^{-1/2}` is the
(non-reduced) Planck mass. This `decay_constant_ratio = f_a / M_Pl` is a
required input of the spindown driver.
"""

"""
    EvolutionDiagnostics

Snapshot recorded at every right-hand-side evaluation: time, amplitudes,
complex mode frequencies, reference-unit charges and losses, and (with
spindown) the background state and metric work.
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
end

"""
    EvolutionResult

Output of a driver: the `OrdinaryDiffEq` solution object (state `[A₁, A₂]` or
`[A₁, A₂, b, j]`), the diagnostics log, and the final radial solution.
"""
struct EvolutionResult{S, T<:Real, R<:RadialSolution}
    ode_solution::S
    diagnostics::Vector{EvolutionDiagnostics{T}}
    final_radial_solution::R
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
end

"""Print one progress line per right-hand-side evaluation when the cache is verbose."""
function _report_step(cache::EvolutionCache, time, amplitudes, solution::RadialSolution, state_text::String)
    cache.verbose || return nothing
    println("t = ", time, "  A = ", amplitudes, "  Im ω = ", growth_rates(solution), "  newton = ", solution.iterations, state_text)
    flush(stdout)
    return nothing
end

function _check_two_modes(problem::RadialProblem)
    n_modes(problem) == 2 || throw(ArgumentError("amplitude evolution needs exactly two bound modes (energy and angular momentum give two equations); got $(n_modes(problem))"))
end

"""Radial solve warm-started from the cache's previous solution."""
function _solve_cached!(cache::EvolutionCache, problem::RadialProblem, amplitudes)
    solution = solve_radial(problem, amplitudes; guess=cache.last_solution, tol=cache.newton_tol)
    cache.last_solution = solution
    return solution
end

# --- fixed background -----------------------------------------------------------

"""
    amplitude_rates(solution)

`dA/dt` on a fixed background from `[∂_A E; ∂_A J] A' = −[L_E; L_J]`.
"""
function amplitude_rates(solution::RadialSolution)
    _check_two_modes(solution.problem)
    return -(charge_jacobian(solution) \ collect(boundary_losses(solution)))
end

function _fixed_background_rhs!(rates, amplitudes, cache::EvolutionCache, time)
    solution = _solve_cached!(cache, cache.template, amplitudes)
    rates .= amplitude_rates(solution)
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(amplitudes), mode_frequencies(solution),
                                                  cloud_charges(solution), boundary_losses(solution),
                                                  1.0, cache.template.background.a, 0.0))
    _report_step(cache, time, amplitudes, solution, "")
    return nothing
end

"""
    evolve_fixed_background(problem, amplitudes, time_span; newton_tol=1e-10, verbose=false, kwargs...)

Integrate the two amplitudes with Tsit5 on a frozen Kerr background.
`verbose` prints one line per right-hand-side evaluation. Keyword arguments
other than these are passed to `solve` (e.g. `reltol`, `abstol`, `saveat`, `dtmax`).
"""
function evolve_fixed_background(problem::RadialProblem, amplitudes::AbstractVector, time_span;
                                 newton_tol=1e-10, verbose=false, reltol=1e-6, abstol=1e-12, kwargs...)
    _check_two_modes(problem)
    cache = EvolutionCache(problem, problem.mass_coupling, 0.0, false, nothing, EvolutionDiagnostics{Float64}[], float(newton_tol), verbose)
    ode = ODEProblem(_fixed_background_rhs!, collect(float.(amplitudes)), time_span, cache)
    ode_solution = solve(ode, Tsit5(); reltol, abstol, kwargs...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution)
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
    metric_work(solution, state, state_rate)

Rate of work done on the cloud by the time-dependent metric, in reference
units per unit `τ`,

    w = ½ ∫ [ (∂_τ g^{αβ}) ∂_αΦ ∂_βΦ + (½ g^{αβ}∂_αΦ∂_βΦ + V) ∂_τ ln(−g) ] √−g d³x,

with `∂_τ` taken at fixed Boyer–Lindquist coordinates in reference units along
the direction `(b', j')` given by `state_rate`.
"""
function metric_work(solution::RadialSolution, state::BackgroundState, state_rate::BackgroundState)
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
    background_charge_derivatives(cache, solution, state; relative_step=1e-5)

Total derivatives `d Q_ref / d b` and `d Q_ref / d j` of the reference-unit
charges at fixed amplitudes and fixed tortoise mesh (in units of the current
mass), returned as a pair of `Charges`. By the implicit-function theorem,

    dQ/db = ∂Q/∂b|_u + (∂Q/∂u) du/db,    du/db = −J⁻¹ ∂F/∂b|_u,

where `∂F/∂b|_u` is exact (ForwardDiff through the background, mass coupling
and mesh offsets: the outer-boundary rows vary as `e^{−κr}` with
`d(κr)/dα ≈ 10³`, far too steeply for finite differences), `∂Q/∂b|_u` is a
central difference of the charges over the background (benign), `J` is the
factorised Jacobian of the converged solution, and the last term is one dual
pass of the charges along `du/db`, exactly as in `charge_jacobian`.
"""
function background_charge_derivatives(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                                       relative_step=1e-6)
    problem, layout, unknowns, scale, amplitudes = solution.problem, solution.layout, solution.unknowns, solution.scale,
                                                   solution.amplitudes
    function derivative_along(direction::BackgroundState, step)
        plus = BackgroundState(state.mass + step * direction.mass, state.angular_momentum + step * direction.angular_momentum)
        minus = BackgroundState(state.mass - step * direction.mass, state.angular_momentum - step * direction.angular_momentum)
        problem_plus, problem_minus = current_problem(cache, plus), current_problem(cache, minus)
        residual_partial = ForwardDiff.derivative(0.0) do epsilon
            perturbed = BackgroundState(state.mass + epsilon * direction.mass,
                                        state.angular_momentum + epsilon * direction.angular_momentum)
            residual(current_problem(cache, perturbed), layout, unknowns, amplitudes, scale)
        end
        charges_partial = (collect(reference_charges(cloud_charges(problem_plus, layout, unknowns, scale), plus)) .-
                           collect(reference_charges(cloud_charges(problem_minus, layout, unknowns, scale), minus))) ./ (2step)
        unknown_direction = -(solution.jacobian_factorization \ residual_partial)
        charges_via_unknowns = ForwardDiff.derivative(0.0) do epsilon
            collect(reference_charges(cloud_charges(problem, layout, unknowns .+ epsilon .* unknown_direction, scale), state))
        end
        total = charges_partial .+ charges_via_unknowns
        return Charges(total[1], total[2])
    end
    mass_step = relative_step * state.mass
    spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
    wrt_mass = derivative_along(BackgroundState(1.0, 0.0), mass_step)
    wrt_spin = derivative_along(BackgroundState(0.0, 1.0), spin_step)
    return wrt_mass, wrt_spin
end

"""
    BalanceTerms

Everything entering the evolving-background balance at one instant, in
reference units: amplitude derivatives of the charges, background derivatives
`∂_b`, `∂_j` of the charges at fixed amplitudes, losses, metric work, and the
background rates driven by the horizon fluxes.
"""
struct BalanceTerms{T<:Real}
    charge_jacobian::Matrix{T}
    charges_wrt_mass::Charges{T}
    charges_wrt_angular_momentum::Charges{T}
    losses::Charges{T}
    metric_work::T
    state_rate::BackgroundState{T}
end

"""
    balance_terms(cache, solution, state; resolve_background=false, relative_step=1e-4)

Assemble the balance terms. The background rates are the horizon fluxes times
`(f_a/M_Pl)²` (see the module docstring). Background derivatives of the charges
(amplitudes and tortoise mesh in units of the current mass held fixed) come
from `background_charge_derivatives`; with `resolve_background=true` they are
instead central finite differences of the charges over four re-solved radial
problems (the slow reference used for validation). The metric work is included
only when the cache says so (see `evolve_with_spindown`); otherwise it is zero
and the balance states that the cloud's energy changes only through its
boundary fluxes, the leading-order model of the literature.
"""
function balance_terms(cache::EvolutionCache, solution::RadialSolution, state::BackgroundState;
                       resolve_background=false, relative_step=1e-4)
    amplitudes = solution.amplitudes
    inner_flux = outward_flux(solution, 1)
    state_rate = BackgroundState(-cache.coupling_to_hole * inner_flux.energy,
                                 -cache.coupling_to_hole * state.mass * inner_flux.angular_momentum)

    if resolve_background
        function charges_at(perturbed::BackgroundState)
            perturbed_solution = solve_radial(current_problem(cache, perturbed), amplitudes; guess=solution, tol=cache.newton_tol)
            reference_charges(cloud_charges(perturbed_solution), perturbed)
        end
        mass_step = relative_step * state.mass
        spin_step = relative_step * max(abs(state.angular_momentum), state.mass^2)
        wrt_mass = (1 / 2mass_step) * (charges_at(BackgroundState(state.mass + mass_step, state.angular_momentum)) -
                                       charges_at(BackgroundState(state.mass - mass_step, state.angular_momentum)))
        wrt_spin = (1 / 2spin_step) * (charges_at(BackgroundState(state.mass, state.angular_momentum + spin_step)) -
                                       charges_at(BackgroundState(state.mass, state.angular_momentum - spin_step)))
    else
        wrt_mass, wrt_spin = background_charge_derivatives(cache, solution, state)
    end

    jacobian = charge_jacobian(solution) .* [state.mass, state.mass^2]
    work = cache.include_metric_work ? metric_work(solution, state, state_rate) : 0.0
    BalanceTerms(jacobian, wrt_mass, wrt_spin, reference_losses(boundary_losses(solution), state), work, state_rate)
end

"""
    amplitude_rates(terms::BalanceTerms)

`A'` from `[∂_A E; ∂_A J] A' = −L + (w, 0) − (∂_b Q b' + ∂_j Q j')`.
"""
function amplitude_rates(terms::BalanceTerms)
    background_drift = terms.state_rate.mass * terms.charges_wrt_mass +
                       terms.state_rate.angular_momentum * terms.charges_wrt_angular_momentum
    right_hand_side = -collect(terms.losses) .+ [terms.metric_work, 0.0] .- collect(background_drift)
    return terms.charge_jacobian \ right_hand_side
end

function _spindown_rhs!(rates, state_vector, cache::EvolutionCache, time)
    amplitudes = state_vector[1:2]
    state = BackgroundState(state_vector[3], state_vector[4])
    solution = _solve_cached!(cache, current_problem(cache, state), amplitudes)
    terms = balance_terms(cache, solution, state)
    rates[1:2] = amplitude_rates(terms)
    rates[3] = terms.state_rate.mass
    rates[4] = terms.state_rate.angular_momentum
    push!(cache.diagnostics, EvolutionDiagnostics(time, copy(amplitudes), mode_frequencies(solution),
                                                  reference_charges(cloud_charges(solution), state), terms.losses,
                                                  state.mass, dimensionless_spin(state), terms.metric_work))
    _report_step(cache, time, amplitudes, solution, "  b = $(state.mass)  j = $(state.angular_momentum)  w = $(terms.metric_work)")
    return nothing
end

"""
    evolve_with_spindown(problem, amplitudes, time_span, decay_constant_ratio;
                         include_metric_work=false, newton_tol=1e-10, verbose=false, kwargs...)

Integrate `[A₁, A₂, b, j]` with Tsit5, starting from `b = 1` and
`j = a` of `problem.background`. `problem` defines the reference units: its
`mass_coupling` is `α₀ = M₀ μ` and its mesh is reused in units of the current
mass at every step. `decay_constant_ratio = f_a / M_Pl` (non-reduced Planck
mass) sets how strongly the cloud's horizon fluxes change the hole.

`include_metric_work` adds the `metric_work` term to the energy balance. It is
off by default: in its present Boyer–Lindquist-fixed form the term grows like
`1/(r − r₊)` towards the horizon and so depends on the inner endpoint (see
NOTES.md); it is also of relative order `(f_a/M_Pl)² E`, like the background
drift terms, so the default is the leading-order balance.
"""
function evolve_with_spindown(problem::RadialProblem, amplitudes::AbstractVector, time_span, decay_constant_ratio::Real;
                              include_metric_work=false, newton_tol=1e-10, verbose=false, reltol=1e-6, abstol=1e-12, kwargs...)
    _check_two_modes(problem)
    cache = EvolutionCache(problem, problem.mass_coupling, float(decay_constant_ratio)^2, include_metric_work, nothing,
                           EvolutionDiagnostics{Float64}[], float(newton_tol), verbose)
    initial = vcat(collect(float.(amplitudes)), [1.0, problem.background.a])
    ode = ODEProblem(_spindown_rhs!, initial, time_span, cache)
    ode_solution = solve(ode, Tsit5(); reltol, abstol, kwargs...)
    EvolutionResult(ode_solution, cache.diagnostics, cache.last_solution)
end
