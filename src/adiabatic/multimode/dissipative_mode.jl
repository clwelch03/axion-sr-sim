"""
    dissipative_source(primary_values, secondary_values, radius, spin, alpha,
                       dissipative_modes, polar_values, azimuthal_values, grid)

Return the source of the m = 0 field at 2ω₁ - ω₂ for each mode at one radius,
as a dictionary keyed by mode: the projection of
`alpha^2 * (r^2 + a^2 cos^2θ) * [1 - cos ϕ₁]₂ * conj(P₂)`, where `[1 - cos ϕ₁]₂`
is [`mass_deficit_harmonic`](@ref) with `azimuthal_number = 2` and
`P₂ = Σ_l R_l Y_l2` is the secondary's polar profile.

The source enters the radial equation in the slot of the nonlinear source `N`
(Omiya et al. 2211.01949 Eqs. (51)-(52), with the factor mu^2).
"""
function dissipative_source(primary_values, secondary_values, radius, spin, alpha,
                            dissipative_modes, polar_values, azimuthal_values, grid)
    # rho^2 = r^2 + a^2 cos^2θ at each polar node.
    rho_squared = radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2

    harmonic = mass_deficit_harmonic(
        primary_values, polar_values, azimuthal_values, grid, 2,
    )
    secondary_profile = sum(
        value .* polar_values[mode] for (mode, value) in secondary_values
    )
    weighted_term = alpha^2 .* grid.polar_weights .* rho_squared .*
                    harmonic .* conj.(secondary_profile)

    return Dict(mode => sum(weighted_term .* polar_values[mode])
                for mode in dissipative_modes)
end

"""
    dissipative_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                              spin, alpha, secondary_frequency, dissipative_frequency,
                              secondary_modes, dissipative_modes, couplings,
                              polar_values, azimuthal_values, grid; include_source=true)

Write the radial derivative of the secondary and the m = 0 field together into
`derivative` and return it.

The state is `[R₂; Q₂; R₀; Q₀]`: the secondary in the layout of
[`secondary_tortoise_rhs!`](@ref), followed by the m = 0 field in
`dissipative_modes` order (all with m = 0, frequency `dissipative_frequency`,
real). Integrating both together makes the secondary available at every RK4
stage. The m = 0 field sees the primary's averaged background, as the
secondary does, plus [`dissipative_source`](@ref) unless `include_source=false`.
"""
function dissipative_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                                   spin, alpha, secondary_frequency, dissipative_frequency,
                                   secondary_modes, dissipative_modes, couplings,
                                   polar_values, azimuthal_values, grid; include_source=true)
    all(mode -> last(mode) == 0, dissipative_modes) ||
        throw(ArgumentError("dissipative modes must have m = 0"))
    secondary_count = length(secondary_modes)
    count = length(dissipative_modes)
    offset = 2secondary_count
    length(state) == length(derivative) == offset + 2count ||
        throw(DimensionMismatch("expected [R2; Q2; R0; Q0]"))

    # The secondary block is unchanged by the m = 0 field.
    secondary_tortoise_rhs!(
        view(derivative, 1:offset), view(state, 1:offset), tortoise_radius,
        primary_values, spin, alpha, secondary_frequency, secondary_modes, couplings,
        polar_values, azimuthal_values, grid,
    )

    geometry = tortoise_geometry(tortoise_radius, spin)
    values = state[offset+1:offset+count]
    slopes = state[offset+count+1:offset+2count]

    # Averaged background alpha^2 V R0, as for the secondary.
    mass_deficit_average = averaged_mass_deficit(
        primary_values, polar_values, azimuthal_values, grid,
    )
    background = alpha^2 .* background_coupling_matrix(
        mass_deficit_average, geometry.radius, spin, dissipative_modes, polar_values, grid,
    ) * values

    source = Dict(mode => zero(eltype(state)) for mode in dissipative_modes)
    if include_source
        secondary_values = Dict(
            mode => state[index] for (index, mode) in enumerate(secondary_modes)
        )
        source = dissipative_source(
            primary_values, secondary_values, geometry.radius, spin, alpha,
            dissipative_modes, polar_values, azimuthal_values, grid,
        )
    end

    # linear_tortoise_rhs! requires m > 0, so assemble the m = 0 operator here.
    for (index, (l, m)) in enumerate(dissipative_modes)
        coefficients = tortoise_operator_coefficients(
            geometry, spin, alpha, dissipative_frequency, l, m,
        )
        coupled_value = sum(
            couplings[(l, other_l, m)] * values[other_index]
            for (other_index, (other_l, _)) in enumerate(dissipative_modes)
        )
        derivative[offset+index] = slopes[index]
        derivative[offset+count+index] =
            -coefficients.slope_coefficient * slopes[index] -
            coefficients.local_coefficient * values[index] -
            coefficients.coupling_scale * coupled_value -
            coefficients.source_scale * (background[index] + source[(l, m)])
    end

    return derivative
end

"""
    dissipative_boundary_state(side, position, spin, alpha, frequency, modes,
                               couplings, amplitudes; horizon_order=3, outer_order=7)

Return the m = 0 field's boundary state `[R₀; Q₀]` in `modes` order, with
`Q = dR/dr_*`, at tortoise `position` on `side = :horizon` or `:outer`.

Ingoing at the horizon; at the outer boundary decaying when
`frequency < alpha`, outgoing otherwise. [`boundary_state`](@ref) requires
m > 0, so this calls the series helpers directly with m = 0 and the actual
`frequency` 2ω₁ - ω₂. Supply one amplitude per mode.
"""
function dissipative_boundary_state(side, position, spin, alpha, frequency, modes,
                                    couplings, amplitudes; horizon_order=3, outer_order=7)
    all(mode -> last(mode) == 0, modes) ||
        throw(ArgumentError("dissipative modes must have m = 0"))
    l_values = first.(modes)

    if side === :horizon
        data = horizon_boundary_data(
            position, spin, alpha, frequency, l_values, 0, couplings, amplitudes;
            order=horizon_order,
        )
    elseif side === :outer
        kind = real(frequency) < alpha ? :bound : :outgoing
        data = outer_boundary_data(
            position, spin, alpha, frequency, l_values, 0, couplings, amplitudes;
            kind, order=outer_order,
        )
    else
        throw(ArgumentError("side must be :horizon or :outer"))
    end

    return vcat(data.values, data.tortoise_slopes)
end

"""
    dissipative_matching_residual(horizon_amplitudes, outer_amplitudes, frequency,
                                  secondary_solution, primary_profiles, primary_modes,
                                  spin, alpha, secondary_modes, dissipative_modes,
                                  couplings, polar_values, azimuthal_values, grid;
                                  include_source=true, horizon_order=3, outer_order=7,
                                  save_profile=false)

Return the mismatch `[R_horizon - R_outer; Q_horizon - Q_outer]` of the m = 0
field at the primary's matching point.

Each side starts the secondary from `secondary_solution` (from
[`solve_secondary_mode`](@ref)) and the m = 0 field from
[`dissipative_boundary_state`](@ref) with the given amplitudes, then integrates
both together with [`dissipative_tortoise_rhs!`](@ref) on the primary's domain
at twice its step. Only the m = 0 part enters the residual. With
`save_profile=true`, also return both integrated branches.
"""
function dissipative_matching_residual(
    horizon_amplitudes, outer_amplitudes, frequency, secondary_solution,
    primary_profiles, primary_modes, spin, alpha, secondary_modes, dissipative_modes,
    couplings, polar_values, azimuthal_values, grid;
    include_source=true, horizon_order=3, outer_order=7, save_profile=false,
)
    primary_from_horizon = primary_profiles.from_horizon
    primary_from_outer = primary_profiles.from_outer
    inner = first(primary_from_horizon.positions)
    match = last(primary_from_horizon.positions)
    outer = first(primary_from_outer.positions)

    secondary_m = last(first(secondary_modes))
    secondary_frequency = secondary_solution.frequency
    initial_state(side, position, secondary_amplitudes, amplitudes) = vcat(
        boundary_state(
            side, position, spin, alpha, secondary_frequency, secondary_modes, couplings,
            secondary_amplitudes; fundamental_m=secondary_m, horizon_order, outer_order,
        ),
        dissipative_boundary_state(
            side, position, spin, alpha, frequency, dissipative_modes, couplings,
            amplitudes; horizon_order, outer_order,
        ),
    )

    # Right-hand side that reads the primary from one saved branch.
    rhs_reading(primary_branch) = (derivative, state, position) ->
        dissipative_tortoise_rhs!(
            derivative, state, position,
            primary_values_at(primary_branch, position, primary_modes),
            spin, alpha, secondary_frequency, frequency, secondary_modes,
            dissipative_modes, couplings, polar_values, azimuthal_values, grid;
            include_source,
        )

    # Twice the primary spacing puts each RK4 half-step on a primary sample,
    # as in secondary_matching_residual.
    primary_spacing = abs(primary_from_horizon.positions[2] - primary_from_horizon.positions[1])
    step = 2 * primary_spacing * (1 + 1e-9)

    from_horizon = integrate_radial(
        rhs_reading(primary_from_horizon),
        initial_state(:horizon, inner, secondary_solution.horizon_amplitudes, horizon_amplitudes),
        inner, match; step, save_profile,
    )
    from_outer = integrate_radial(
        rhs_reading(primary_from_outer),
        initial_state(:outer, outer, secondary_solution.outer_amplitudes, outer_amplitudes),
        outer, match; step, save_profile,
    )

    # The m = 0 field follows the secondary's 2N entries.
    dissipative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2length(dissipative_modes)
    if save_profile
        residual = from_horizon.states[dissipative_part, end] - from_outer.states[dissipative_part, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon[dissipative_part] - from_outer[dissipative_part]
end

"""
    solve_dissipative_mode(frequency, secondary_solution, primary_profiles,
                           primary_modes, spin, alpha, secondary_modes,
                           dissipative_modes, couplings, polar_values,
                           azimuthal_values, grid; series_options...)

Solve for the m = 0 field at `frequency` (2 Re ω₁ - Re ω₂) driven by the
primary and the secondary, and return its amplitudes, profiles, and fluxes.

The residual of [`dissipative_matching_residual`](@ref) is affine in the 2N
boundary amplitudes `c`: `residual(c) = r0 + M * c`. Build `r0` with the source
and all amplitudes zero, and column i of `M` with the source off and unit
amplitude i, then solve `M * c = -r0`. Fluxes are in units `Fa² M` per unit
(2,2) outer amplitude of the secondary squared; positive means energy lost.
`condition_number` is that of `M` with columns scaled to unit norm.
"""
function solve_dissipative_mode(frequency, secondary_solution, primary_profiles,
                                primary_modes, spin, alpha, secondary_modes,
                                dissipative_modes, couplings, polar_values,
                                azimuthal_values, grid; series_options...)
    count = length(dissipative_modes)
    residual_at(amplitudes; include_source) = dissipative_matching_residual(
        amplitudes[1:count], amplitudes[count+1:end], frequency, secondary_solution,
        primary_profiles, primary_modes, spin, alpha, secondary_modes, dissipative_modes,
        couplings, polar_values, azimuthal_values, grid; include_source, series_options...,
    )

    source_residual = residual_at(zeros(ComplexF64, 2count); include_source=true)
    matching_matrix = hcat([
        residual_at([i == j ? 1.0 + 0im : 0im for j in 1:2count]; include_source=false)
        for i in 1:2count
    ]...)
    amplitudes = -(matching_matrix \ source_residual)

    # A unit outer amplitude grows by many orders of magnitude on its way in to
    # the matching point, so the raw columns differ wildly in size. Rescaling a
    # column only rescales its unknown; the condition number of the rescaled
    # matrix measures how close the problem is to singular (a resonance).
    column_norms = [norm(column) for column in eachcol(matching_matrix)]
    condition_number = cond(matching_matrix ./ transpose(column_norms))

    profiles = dissipative_matching_residual(
        amplitudes[1:count], amplitudes[count+1:end], frequency, secondary_solution,
        primary_profiles, primary_modes, spin, alpha, secondary_modes, dissipative_modes,
        couplings, polar_values, azimuthal_values, grid; save_profile=true, series_options...,
    )

    # Loss into the horizon is minus the outward current at the inner endpoint.
    dissipative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2count
    horizon_flux = -outward_secondary_flux(
        profiles.from_horizon.states[dissipative_part, 1],
        first(profiles.from_horizon.positions), frequency, spin, count,
    )
    outer_flux = outward_secondary_flux(
        profiles.from_outer.states[dissipative_part, 1],
        first(profiles.from_outer.positions), frequency, spin, count,
    )

    return (;
        horizon_amplitudes=amplitudes[1:count], outer_amplitudes=amplitudes[count+1:end],
        residual_norm=norm(profiles.residual), source_residual_norm=norm(source_residual),
        condition_number, profiles, horizon_flux, outer_flux,
    )
end
