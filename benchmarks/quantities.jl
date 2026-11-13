# Observables compared by the regression check. Included at the end of the
# `BenchmarkHarness` module in common.jl.
#
# A `Quantity` is one named group of numbers (a norm, or a few fixed-index entries,
# or all entries of a small vector). It carries a tolerance `class` (looked up in
# tolerances.toml) and an `anchor`, an absolute size below which differences are
# not significant. The comparison measure is
#
#     relative difference = ‖new − reference‖₂ / max(‖reference‖₂, anchor)
#
# `anchor = 0` for quantities that do not suffer cancellation. For a quantity that
# is a difference of nearly equal terms (losses, rates, work), the anchor is the
# size of the terms, so that the tolerance is relative to the inputs and not to the
# nearly cancelled result. The anchor is computed when the reference is written and
# stored with it.

struct Quantity
    values::Vector{Float64}
    class::String
    anchor::Float64
end

Quantity(values, class::String; anchor=0.0) = Quantity(collect(Float64.(values)), class, Float64(anchor))

const QuantityTable = Vector{Pair{String, Quantity}}

"""Number of sampled entries of a long vector kept in the reference."""
const SAMPLE_COUNT = 12

add!(table::QuantityTable, name::String, quantity::Quantity) = push!(table, name => quantity)

"""Complex numbers as interleaved real and imaginary parts."""
interleave(values::AbstractVector{<:Complex}) = collect(Iterators.flatten((real(value), imag(value)) for value in values))

# --- groups ----------------------------------------------------------------------

"""Residual and Jacobian at the converged and at the displaced unknowns."""
function residual_and_jacobian_quantities!(table::QuantityTable, state::BenchmarkState)
    solution = state.solution
    displaced = displaced_unknowns(solution)
    indices = sample_indices(length(displaced), SAMPLE_COUNT)
    displaced_residual = residual_at(state, displaced)
    converged_residual = residual_at(state, solution.unknowns)
    # at the root every row is roundoff-sized; compare it to the size of the displaced residual
    add!(table, "residual.converged_max_abs", Quantity([norm(converged_residual, Inf)], "residual_converged";
                                                       anchor=norm(displaced_residual, Inf)))
    add!(table, "residual.displaced_norm", Quantity([norm(displaced_residual)], "residual"))
    add!(table, "residual.displaced_entries", Quantity(displaced_residual[indices], "residual"))
    matrix = jacobian_at(state, displaced)
    vector = fixed_vector(length(displaced))
    product = matrix * vector
    add!(table, "jacobian.stored_entries", Quantity([nnz(matrix)], "exact"))
    add!(table, "jacobian.frobenius_norm", Quantity([norm(matrix)], "jacobian"))
    add!(table, "jacobian.times_vector_norm", Quantity([norm(product)], "jacobian"))
    add!(table, "jacobian.times_vector_entries", Quantity(product[indices], "jacobian"; anchor=1e-3 * norm(product, Inf)))
end

"""Converged `solve_radial` after the fixed perturbation of amplitudes and background."""
function perturbed_solve_quantities!(table::QuantityTable, state::BenchmarkState; newton_tolerance, guess)
    solution = perturbed_solve(state; tolerance=newton_tolerance, guess)
    solution.converged || error("perturbed solve did not converge on state '$(state.spec.name)'")
    indices = sample_indices(length(solution.unknowns), SAMPLE_COUNT)
    charges = cloud_charges(solution)
    add!(table, "solve.mode_frequencies", Quantity(interleave(mode_frequencies(solution)), "solve_frequency"))
    add!(table, "solve.unknown_entries", Quantity(solution.unknowns[indices], "solve_unknown";
                                                  anchor=1e-3 * norm(solution.unknowns, Inf)))
    add!(table, "solve.cloud_charges", Quantity([charges.energy, charges.angular_momentum], "solve_charge"))
    add!(table, "solve.iterations", Quantity([solution.iterations], "reported"))
    return solution
end

"""Size of the larger of the inner and outer fluxes: the terms whose difference is a loss."""
function edge_flux_scales(solution::RadialSolution)
    last_node = solution.layout.n_nodes
    inner, outer = outward_flux(solution, 1), outward_flux(solution, last_node)
    energy = max(abs(inner.energy), abs(outer.energy))
    angular_momentum = max(abs(inner.angular_momentum), abs(outer.angular_momentum))
    action = max(maximum(abs, action_flux(solution, 1)), maximum(abs, action_flux(solution, last_node)))
    return (charges=hypot(energy, angular_momentum), action=action)
end

"""Charges, actions and losses at the state itself."""
function charge_quantities!(table::QuantityTable, state::BenchmarkState)
    solution = state.solution
    charges, losses = cloud_charges(solution), boundary_losses(solution)
    scales = edge_flux_scales(solution)
    add!(table, "charges.cloud_charges", Quantity([charges.energy, charges.angular_momentum], "charge"))
    add!(table, "charges.phase_actions", Quantity(phase_actions(solution), "charge"))
    # losses = outer flux − inner flux: the anchor is the size of the two fluxes
    add!(table, "charges.boundary_losses", Quantity(collect(losses), "loss"; anchor=scales.charges))
    add!(table, "charges.action_losses", Quantity(action_losses(solution), "loss"; anchor=scales.action))
end

"""Amplitude response, charge and action Jacobians."""
function response_quantities!(table::QuantityTable, state::BenchmarkState)
    solution = state.solution
    directions = amplitude_response(solution)
    indices = sample_indices(size(directions, 1), SAMPLE_COUNT)
    add!(table, "response.direction_norms", Quantity([norm(column) for column in eachcol(directions)], "response"))
    for (mode, column) in enumerate(eachcol(directions))
        add!(table, "response.direction_$(mode)_entries", Quantity(column[indices], "response"; anchor=1e-3 * norm(column, Inf)))
    end
    add!(table, "response.charge_jacobian", Quantity(vec(charge_jacobian(solution; directions)), "response"))
    add!(table, "response.action_jacobian", Quantity(vec(action_jacobian(solution; directions)), "response"))
end

"""Every field of a `BackgroundDerivatives`, with the group name as prefix."""
function derivative_quantities!(table::QuantityTable, prefix::String, derivatives::BackgroundDerivatives)
    for field in fieldnames(BackgroundDerivatives)
        value = getfield(derivatives, field)
        values = value isa Charges ? [value.energy, value.angular_momentum] : value
        add!(table, "$prefix.$field", Quantity(values, "derivative"))
    end
end

"""`background_charge_derivatives`, `balance_terms` and the rates they give."""
function balance_quantities!(table::QuantityTable, state::BenchmarkState; newton_tolerance)
    solution, mass = state.solution, state.background.mass
    cache = evolution_cache(state; newton_tolerance)
    derivative_quantities!(table, "background_derivatives", background_charge_derivatives(cache, solution, state.background))
    terms = balance_terms(cache, solution, state.background)
    scales = edge_flux_scales(solution)
    derivatives, rate = terms.derivatives, terms.state_rate
    add!(table, "balance.charge_jacobian", Quantity(vec(terms.charge_jacobian), "response"))
    add!(table, "balance.action_jacobian", Quantity(vec(terms.action_jacobian), "response"))
    derivative_quantities!(table, "balance.derivatives", derivatives)
    # losses and work are differences of nearly equal terms; anchors are the sizes of those terms (reference units)
    add!(table, "balance.losses", Quantity(collect(terms.losses), "loss"; anchor=scales.charges))
    add!(table, "balance.action_losses", Quantity(terms.action_losses, "loss"; anchor=mass * scales.action))
    work_terms = abs(derivatives.charges_wrt_mass.energy * rate.mass) +
                 abs(derivatives.charges_wrt_angular_momentum.energy * rate.angular_momentum) +
                 sum(abs.(real.(mode_frequencies(solution)) .* derivatives.actions_wrt_mass)) * abs(rate.mass) +
                 sum(abs.(real.(mode_frequencies(solution)) .* derivatives.actions_wrt_angular_momentum)) * abs(rate.angular_momentum)
    add!(table, "balance.metric_work", Quantity([terms.metric_work], "loss"; anchor=work_terms))
    add!(table, "balance.state_rate", Quantity([rate.mass, rate.angular_momentum], "rate_background"))
    # amplitude rates solve [∂_A N] A' = −L_N − drift; anchor on |[∂_A N]⁻¹| applied to the sizes of the two terms
    drift_size = abs(rate.mass) .* abs.(derivatives.actions_wrt_mass) .+
                 abs(rate.angular_momentum) .* abs.(derivatives.actions_wrt_angular_momentum)
    rate_anchor(jacobian, loss_size, drift) = norm(abs.(inv(jacobian)) * (fill(loss_size, size(jacobian, 1)) .+ drift))
    add!(table, "rates.amplitude_rates", Quantity(amplitude_rates(terms, cache), "rate";
                                                  anchor=rate_anchor(terms.action_jacobian, mass * scales.action, drift_size)))
    add!(table, "rates.fixed_background_amplitude_rates", Quantity(amplitude_rates(solution), "rate";
                                                                   anchor=rate_anchor(action_jacobian(solution), scales.action, 0)))
end

"""Short `evolve_with_spindown` segment: trajectory, step and right-hand-side counts."""
function evolution_quantities!(table::QuantityTable, state::BenchmarkState; time_span=(1.5e9, 1.9e9), newton_tolerance)
    result, rhs_calls = evolution_segment(state; time_span, newton_tolerance)
    trajectory = amplitude_trajectory(result)
    add!(table, "evolution.accepted_steps", Quantity([length(trajectory.times) - 1], "reported"))
    add!(table, "evolution.rhs_calls", Quantity([rhs_calls], "reported"))
    add!(table, "evolution.rhs_diagnostics_logged", Quantity([length(result.diagnostics)], "reported"))
    add!(table, "evolution.times", Quantity(trajectory.times, "trajectory_time"))
    for (mode, amplitudes) in enumerate(eachcol(trajectory.amplitudes))
        add!(table, "evolution.amplitude_$mode", Quantity(amplitudes, "trajectory_amplitude"))
    end
    add!(table, "evolution.black_hole_mass", Quantity(trajectory.black_hole_mass, "trajectory_background"))
    add!(table, "evolution.black_hole_angular_momentum", Quantity(trajectory.black_hole_angular_momentum, "trajectory_background"))
    return result
end

"""
    collect_quantities(state; newton_tolerance=state.spec.run_tolerance, evolution=state.spec.include_evolution,
                       perturbed_guess=state.solution)

All compared observables of one state as a `QuantityTable` with names `"<state>/<group>.<item>"`.
`newton_tolerance` is used by the perturbed solve and inside the evolution cache.
"""
function collect_quantities(state::BenchmarkState; newton_tolerance=state.spec.run_tolerance, evolution=state.spec.include_evolution,
                            perturbed_guess=state.solution)
    table = QuantityTable()
    residual_and_jacobian_quantities!(table, state)
    perturbed_solve_quantities!(table, state; newton_tolerance, guess=perturbed_guess)
    charge_quantities!(table, state)
    response_quantities!(table, state)
    balance_quantities!(table, state; newton_tolerance)
    evolution && evolution_quantities!(table, state; newton_tolerance)
    return [state.spec.name * "/" * name => quantity for (name, quantity) in table]
end

# --- comparison ------------------------------------------------------------------

"""Relative difference of `new` from `reference` as defined at the top of this file."""
function relative_difference(new::AbstractVector, reference::AbstractVector, anchor::Real)
    length(new) == length(reference) || return Inf
    denominator = max(norm(reference), anchor)
    difference = norm(new .- reference)
    return denominator > 0 ? difference / denominator : (difference == 0 ? 0.0 : Inf)
end
