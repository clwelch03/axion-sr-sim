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
                             cloud_frequencies, cloud_m, couplings,
                             source_basis; potential=:cosine)

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
Supply a basis built for the same cloud ordering and angular/phase grids,
with mode labels in exactly the state order. Build it once before radial
integration and reuse it at every step. Select :cosine, :quartic, or :sextic
potential.
"""
function nonlinear_two_phase_rhs!(derivative, state, position, spin, alpha,
                                  components, cloud_frequencies, cloud_m,
                                  couplings, source_basis::TwoPhaseSourceBasis; potential=:cosine)
    mode_labels = [(l, component.q) for component in components for l in component.l_values]
    num_radial_modes = length(mode_labels)
    length(state) == length(derivative) == 2num_radial_modes ||
        throw(DimensionMismatch("expected [R; Q] for all retained pairs"))
    source_basis.mode_labels == mode_labels ||
        throw(ArgumentError("source basis mode order must match the radial state"))

    # Pack R only; the nonlinear potential contains no radial slopes Q.
    geometry = tortoise_geometry(position, spin)
    radial_coefficients = state[1:num_radial_modes]
    parameters = vcat(real.(radial_coefficients), imag.(radial_coefficients),
        [geometry.radius, spin, alpha])
    packed_source = nonlinear_two_phase_source(parameters, source_basis; potential)
    source = complex.(packed_source[1:num_radial_modes],
        packed_source[num_radial_modes+1:2num_radial_modes])
    source_scale = geometry.delta / geometry.radial_factor^2

    offset = 0
    for component in components
        indices = (offset + 1):(offset + length(component.l_values))
        slope_indices = num_radial_modes .+ indices
        operator = component_tortoise_operator(position, spin, alpha, component,
            cloud_frequencies, cloud_m, couplings)

        # Linear mixing stays within this q block; N couples the joint field.
        derivative[indices] .= state[slope_indices]
        derivative[slope_indices] .= -operator.slope_coefficient .* state[slope_indices] -
            operator.profile_matrix * state[indices] -
            source_scale .* source[indices]
        offset += length(component.l_values)
    end
    return derivative
end
