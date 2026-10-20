using LinearAlgebra: Diagonal

"""
    component_tortoise_operator(position, spin, alpha, component,
                               cloud_frequencies, cloud_m, couplings)

Assemble the linear radial equation for one Fourier component in M=1 units.

`position` is the tortoise radius. Supply a validated component and the
same cloud ordering used when constructing it. Supply `couplings` from
`build_cosine_squared_couplings`.

Return `(slope_coefficient, profile_matrix)` defining
`R'' = -slope_coefficient*R' - profile_matrix*R`.
(Derivatives are with respect to tortoise radius.)
"""
function component_tortoise_operator(position, spin, alpha, component,
                                     cloud_frequencies, cloud_m, couplings)
    numbers = component_frequency_and_m(component.q, cloud_frequencies, cloud_m)
    geometry = tortoise_geometry(position, spin)

    radial_coefficients = [tortoise_operator_coefficients(
        geometry, spin, alpha, numbers.frequency, l, numbers.m,
    ) for l in component.l_values]

    angular_matrix = [couplings[(l, other_l, abs(numbers.m))]
        for l in component.l_values, other_l in component.l_values]

    profile_matrix = Diagonal([c.local_coefficient for c in radial_coefficients]) +
        radial_coefficients[1].coupling_scale * angular_matrix

    return (; slope_coefficient=radial_coefficients[1].slope_coefficient, profile_matrix)
end

"""
    component_boundary_state(side, position, spin, alpha, component,
                             cloud_frequencies, cloud_m, couplings, leading_coefficients;
                             horizon_order=3, outer_order=7)

Build linear boundary data for one Fourier component in M=1 units.

Return [R; dR/dr_*], with each half ordered by component.l_values.
`position` is a tortoise radius; `leading_coefficients` gives one leading
asymptotic-series coefficient per retained l. The boundary helpers use these
coefficients to calculate the endpoint field values and slopes.

Choose :horizon for ingoing data or :outer for decaying/outgoing data.
The outer branch is bound when Re(frequency) < alpha, otherwise outgoing.
Require a finite frequency with positive real part. Static components and
negative-frequency representatives are deferred; the outer helper rejects
zero wave number at the mass threshold. Zero and signed m are supported.
"""
function component_boundary_state(side, position, spin, alpha, component,
                                  cloud_frequencies, cloud_m, couplings, leading_coefficients;
                                  horizon_order=3, outer_order=7)
    side in (:horizon, :outer) ||
        throw(ArgumentError("side must be :horizon or :outer"))
    numbers = component_frequency_and_m(component.q, cloud_frequencies, cloud_m)
    isfinite(numbers.frequency) && real(numbers.frequency) > 0 ||
        throw(DomainError(numbers.frequency, "require finite Re(frequency) > 0"))

    block_couplings = Dict(
        (l, other_l, numbers.m) => couplings[(l, other_l, abs(numbers.m))]
        for l in component.l_values for other_l in component.l_values
    )

    if side == :horizon
        data = horizon_boundary_data(
            position, spin, alpha, numbers.frequency, component.l_values,
            numbers.m, block_couplings, leading_coefficients; order=horizon_order)
    else
        kind = real(numbers.frequency) < alpha ? :bound : :outgoing
        data = outer_boundary_data(
            position, spin, alpha, numbers.frequency, component.l_values,
            numbers.m, block_couplings, leading_coefficients; kind, order=outer_order)
    end
    return vcat(data.values, data.tortoise_slopes)
end

"""
    nonlinear_two_phase_rhs!(derivative, state, position, spin, alpha, components,
                             cloud_frequencies, cloud_m, couplings, polar_values,
                             grid, phase1_nodes, phase2_nodes; potential=:cosine)

Calculate the radial derivatives of all retained Fourier components.

`position` is the tortoise radius r_*, in M=1 units.
For each retained (l, q), evolve the radial coefficient R and its slope
Q = dR/dr_* according to

    dR/dr_* = Q
    dQ/dr_* = -p*Q - V*R - Delta/(r^2 + spin^2)^2 * N

For n pairs, `state = [R_1, ..., R_n, Q_1, ..., Q_n]`, with Q = dR/dr_*.
Order them by walking through `components` in the supplied list order,
and through each component's `l_values` in its supplied order.
`derivative` uses exactly the same layout.

Pass separate arrays for `state` and `derivative`, each with two entries
per retained (l, q). The function fills and returns `derivative`,
leaving `state` unchanged.

Compute N from the joint field of all components. Supply validated components
with distinct, nonzero q labels and only one member of each q, -q pair.
Use the same cloud ordering, angular tables, and phase grids as for
`nonlinear_two_phase_source`. Select :cosine, :quartic, or :sextic potential.
"""
function nonlinear_two_phase_rhs!(derivative, state, position, spin, alpha,
                                  components, cloud_frequencies, cloud_m,
                                  couplings, polar_values, grid,
                                  phase1_nodes, phase2_nodes; potential=:cosine)
    pairs = [(l, component.q) for component in components for l in component.l_values]
    count = length(pairs)
    length(state) == length(derivative) == 2count ||
        throw(DimensionMismatch("expected [R; Q] for all retained pairs"))

    # Reconstruct the joint field from R only; the potential contains no Q.
    radial_values = Dict(pair => state[index] for (index, pair) in enumerate(pairs))
    geometry = tortoise_geometry(position, spin)
    source = nonlinear_two_phase_source(radial_values, geometry.radius, spin, alpha,
        pairs, cloud_m, polar_values, grid, phase1_nodes, phase2_nodes; potential)
    source_scale = geometry.delta / geometry.radial_factor^2

    offset = 0
    for component in components
        indices = (offset + 1):(offset + length(component.l_values))
        slope_indices = count .+ indices
        operator = component_tortoise_operator(position, spin, alpha, component,
            cloud_frequencies, cloud_m, couplings)

        # Linear mixing stays within this q block; N couples the joint field.
        derivative[indices] .= state[slope_indices]
        derivative[slope_indices] .= -operator.slope_coefficient .* state[slope_indices] -
            operator.profile_matrix * state[indices] -
            source_scale .* [source[pair] for pair in pairs[indices]]
        offset += length(component.l_values)
    end
    return derivative
end

"""
    two_phase_matching_residual(cloud_frequencies, horizon_coefficients,
        outer_coefficients, spin, alpha, components, cloud_m, couplings,
        polar_values, grid, phase1_nodes, phase2_nodes; kwargs...)

Integrate from both boundaries and return their mismatch at `match`.

Supply one leading boundary-series coefficient per retained (l, q), ordered
by the supplied component list and then each component's l_values.
The frequency vector contains one trial frequency per independent cloud.
This function evaluates a trial; it does not fix normalizations or run Newton.

`inner`, `outer`, and `match` are tortoise radii; require inner < match < outer.
`step` is the maximum radial RK4 spacing. Return [R_inner-R_outer; Q_inner-Q_outer].
With `save_profile=true`, return (from_horizon, from_outer, residual), keeping
both radial integration histories and their separate matching-point values.
"""
function two_phase_matching_residual(cloud_frequencies, horizon_coefficients,
        outer_coefficients, spin, alpha, components, cloud_m, couplings,
        polar_values, grid, phase1_nodes, phase2_nodes;
        potential=:cosine, step=0.05, inner=-100.0, outer=100.0, match=5.0,
        horizon_order=3, outer_order=7, save_profile=false)
    inner < match < outer || throw(ArgumentError("match must lie between endpoints"))
    isempty(components) && throw(ArgumentError("require at least one component"))
    counts = [length(component.l_values) for component in components]
    offsets = cumsum(vcat(0, counts))
    length(horizon_coefficients) == length(outer_coefficients) == offsets[end] ||
        throw(DimensionMismatch("require one boundary coefficient per retained pair"))

    function assemble_boundary(side, position, coefficients)
        blocks = [component_boundary_state(side, position, spin, alpha, component,
            cloud_frequencies, cloud_m, couplings,
            coefficients[offsets[j]+1:offsets[j+1]]; horizon_order, outer_order)
            for (j, component) in enumerate(components)]

        # Each block is [R; Q]. Gather ALL R blocks before ALL Q blocks for the RHS.
        profiles = [block[1:counts[j]] for (j, block) in enumerate(blocks)]
        slopes = [block[counts[j]+1:end] for (j, block) in enumerate(blocks)]
        return vcat(profiles..., slopes...)
    end

    horizon_state = assemble_boundary(:horizon, inner, horizon_coefficients)
    outer_state = assemble_boundary(:outer, outer, outer_coefficients)
    rhs! = (derivative, state, position) -> nonlinear_two_phase_rhs!(
        derivative, state, position, spin, alpha, components, cloud_frequencies,
        cloud_m, couplings, polar_values, grid, phase1_nodes, phase2_nodes; potential)

    from_horizon = integrate_radial(rhs!, horizon_state, inner, match; step, save_profile)
    from_outer = integrate_radial(rhs!, outer_state, outer, match; step, save_profile)

    if save_profile
        residual = from_horizon.states[:, end] - from_outer.states[:, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon - from_outer
end

"""
    unpack_two_phase_trial(frequency_and_boundary_guess, cloud_amplitudes,
        fourier_components; coefficient_scale)

Organize a candidate solution into frequencies, horizon boundary coefficients,
and outer boundary coefficients.

The guess contains, in order:
- Two cloud frequencies.
- All horizon boundary coefficients, divided by coefficient_scale.
- The free outer boundary coefficients, divided by coefficient_scale.

Multiply the stored boundary coefficients by coefficient_scale to recover their
physical values. Insert the two prescribed cloud_amplitudes into the outer
coefficient array, at the first retained angular harmonic of each fundamental
component, q=(1,0) and q=(0,1).

For this example with |211>, |322>, and two driven components (2w_1 - w_2) and
(2w_2 - w_1) we have eight boundary coefficients, two of which (the outer parent cloud
coefficients) we prescribe as A_1 and A_2. Thus we would return

frequencies          = [ω1, ω2]
horizon_coefficients = [h1, h2, h3, h4]
outer_coefficients   = [A1, A2, b3, b4]
"""
function unpack_two_phase_trial(frequency_and_boundary_guess, cloud_amplitudes, fourier_components; coefficient_scale)
    length(cloud_amplitudes) == 2 ||
        throw(DimensionMismatch("require two normalizations"))
    all(value -> value isa Real && isfinite(value) && value > 0, cloud_amplitudes) ||
        throw(ArgumentError("require positive finite real normalizations"))
    coefficient_scale isa Real && isfinite(coefficient_scale) && coefficient_scale > 0 ||
        throw(ArgumentError("require a positive finite real coefficient scale"))

    radial_mode_labels = [(l, component.q) for component in fourier_components for l in component.l_values]
    num_radial_modes = length(radial_mode_labels)
    length(frequency_and_boundary_guess) ==  2 * num_radial_modes ||
        throw(DimensionMismatch("require two frequencies and all free coefficients"))
    fixed = [findfirst(pair -> pair[2] == q, radial_mode_labels) for q in ((1, 0), (0, 1))]
    any(isnothing, fixed) && throw(ArgumentError("require both fundamental components"))
    fixed = Int.(fixed)
    free = setdiff(1:num_radial_modes, fixed)

    # Promote both arguments: amplitude differentiation also needs dual-valued
    # boundary storage when the trial state itself contains ordinary numbers.
    value_type = promote_type(eltype(frequency_and_boundary_guess), eltype(cloud_amplitudes), typeof(coefficient_scale))
    frequencies = value_type.(frequency_and_boundary_guess[1:2])
    horizon = coefficient_scale .* value_type.(frequency_and_boundary_guess[3:count+2])
    outer = zeros(value_type, count)
    outer[fixed] .= cloud_amplitudes
    outer[free] .= coefficient_scale .* frequency_and_boundary_guess[count+3:end]
    return frequencies, horizon, outer
end

"""
    solve_two_phase_mode(frequency_guess, normalizations, spin, alpha,
        components, cloud_m, couplings, polar_values, grid, phase1_nodes,
        phase2_nodes; horizon_guess, outer_guess, kwargs...)

Solve matching at two prescribed positive outer coefficients.

`normalizations[a]` fixes the first retained l in fundamental q=(1,0) or (0,1).
Supply full boundary-coefficient guess vectors in the RHS pair order.
The two fixed entries of `outer_guess` are replaced by `normalizations`.

Solve for two complex frequencies, all horizon coefficients, and free outer
coefficients. Use validated components with distinct nonzero q labels,
one per conjugate pair.
    
Scale coefficients and matching residuals by the smaller normalization;
`tolerance` applies to that scaled residual norm. Keep both fundamentals bound
and every component on its initial positive-frequency boundary branch.
Forward integration keywords to `two_phase_matching_residual`.

Return frequencies, physical boundary coefficients, residuals, and iterations.
This numerical prototype does not derive a nonlinear amplitude-evolution law.
"""
function solve_two_phase_mode(frequency_guess, normalizations, spin, alpha,
        components, cloud_m, couplings, polar_values, grid, phase1_nodes, phase2_nodes;
        horizon_guess, outer_guess, potential=:cosine, tolerance=1e-9,
        max_iterations=15, difference_step=1e-6, progress_callback=nothing,
        integration_options...)
    length(frequency_guess) == length(normalizations) == length(cloud_m) == 2 ||
        throw(DimensionMismatch("require two frequencies, normalizations, and cloud m values"))
    all(value -> value isa Real && isfinite(value) && value > 0, normalizations) ||
        throw(ArgumentError("require positive finite real outer normalizations"))

    pairs = [(l, component.q) for component in components for l in component.l_values]
    count = length(pairs)
    length(horizon_guess) == length(outer_guess) == count ||
        throw(DimensionMismatch("require a full boundary guess at each endpoint"))
    fixed = [findfirst(pair -> pair[2] == q, pairs) for q in ((1, 0), (0, 1))]
    any(isnothing, fixed) && throw(ArgumentError("require both fundamental components"))
    fixed = Int.(fixed)
    free = setdiff(1:count, fixed)
    scale = minimum(normalizations)

    # Newton's vector contains frequencies, scaled horizon coefficients, then
    # scaled free outer coefficients. Reinsert the two fixed outer coefficients.
    unpack_trial = state -> unpack_two_phase_trial(
        state, normalizations, components; coefficient_scale=scale)
    initial = complex.(vcat(frequency_guess, horizon_guess ./ scale, outer_guess[free] ./ scale))

    initial_bound = [real(sum(component.q .* frequency_guess)) < alpha for component in components]
    function admissible(state)
        all(isfinite, state) && all(value -> 0 < real(value) < alpha, state[1:2]) || return false
        frequencies = [real(sum(component.q .* state[1:2])) for component in components]
        return all(value -> value > 0 && value != alpha, frequencies) &&
            (frequencies .< alpha) == initial_bound
    end

    residual = state -> two_phase_matching_residual(unpack_trial(state)...,
        spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes, phase2_nodes; potential, save_profile=false, integration_options...) ./ scale

    # A real Jacobian is required because field reconstruction uses conjugates.
    solution = real_newton(residual, initial; admissible, tolerance, max_iterations,
        difference_step, progress_callback)

    frequencies, horizon_coefficients, outer_coefficients = unpack_trial(solution.state)
    physical_residual = scale .* solution.residual
    return (; frequencies, horizon_coefficients, outer_coefficients,
        residual=physical_residual, residual_norm=norm(physical_residual),
        scaled_residual_norm=norm(solution.residual), iterations=solution.iterations)
end
