"""
    radiative_source(primary_values, secondary_values, radius, spin, alpha,
                     radiative_modes, polar_values, azimuthal_values, grid)

Return the source of the m = 3 field at 2ω₂ - ω₁ for each mode at one radius,
as a dictionary keyed by mode: the projection of
`(1/2) * alpha^2 * (r^2 + a^2 cos^2θ) * conj([sin ϕ₁]₁) * P₂^2`, where
`[sin ϕ₁]₁` is [`sine_harmonic`](@ref) with `azimuthal_number = 1` and
`P₂ = Σ_l R_l Y_l2` is the secondary's polar profile.

The source enters the radial equation in the slot of the nonlinear source `N`,
with the sign of [`dissipative_source`](@ref): `N` collects
`alpha^2 * rho^2 * (ϕ - sin ϕ)`, whose second-order term in ϕ₂ is
`(1/2) sin ϕ₁ ϕ₂^2` (Omiya et al. 2211.01949 Eqs. (51)-(52), with the factor mu^2).
"""
function radiative_source(primary_values, secondary_values, radius, spin, alpha,
                          radiative_modes, polar_values, azimuthal_values, grid)
    # rho^2 = r^2 + a^2 cos^2θ at each polar node.
    rho_squared = radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2

    # The part of sin ϕ₁ varying as exp(+im*ψ₁) is the conjugate of the k = 1 harmonic.
    harmonic = conj.(sine_harmonic(primary_values, polar_values, azimuthal_values, grid, 1))
    secondary_profile = sum(
        value .* polar_values[mode] for (mode, value) in secondary_values
    )
    weighted_term = (alpha^2 / 2) .* grid.polar_weights .* rho_squared .*
                    harmonic .* secondary_profile .^ 2

    return Dict(mode => sum(weighted_term .* polar_values[mode])
                for mode in radiative_modes)
end

"""
    radiative_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                            spin, alpha, secondary_frequency, radiative_frequency,
                            secondary_modes, radiative_modes, couplings,
                            polar_values, azimuthal_values, grid; include_source=true)

Write the radial derivative of the combined state `[R₂; Q₂; R₃; Q₃]` into
`derivative` and return it.

The secondary is carried along only so that its values are available at every
RK4 stage (the source needs them at step midpoints, where no saved sample
exists). Its block is computed by [`secondary_tortoise_rhs!`](@ref) alone and
does not depend on the m = 3 field.

The state is `[R₂; Q₂; R₃; Q₃]`: the secondary in the layout of
[`secondary_tortoise_rhs!`](@ref), followed by the m = 3 field in
`radiative_modes` order (all with m = 3, frequency `radiative_frequency`, real).
The m = 3 field sees the primary's averaged background, as the secondary does,
so its source-free part is [`secondary_tortoise_rhs!`](@ref) with m = 3 and
ω₃. [`radiative_source`](@ref) is added unless `include_source=false`.
"""
function radiative_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                                 spin, alpha, secondary_frequency, radiative_frequency,
                                 secondary_modes, radiative_modes, couplings,
                                 polar_values, azimuthal_values, grid; include_source=true)
    all(mode -> last(mode) == 3, radiative_modes) ||
        throw(ArgumentError("radiative modes must have m = 3"))
    count = length(radiative_modes)
    offset = 2length(secondary_modes)
    length(state) == length(derivative) == offset + 2count ||
        throw(DimensionMismatch("expected [R2; Q2; R3; Q3]"))

    # The secondary block is unchanged by the m = 3 field.
    secondary_tortoise_rhs!(
        view(derivative, 1:offset), view(state, 1:offset), tortoise_radius,
        primary_values, spin, alpha, secondary_frequency, secondary_modes, couplings,
        polar_values, azimuthal_values, grid,
    )
    # The m = 3 block: Kerr operator at ω₃ plus the averaged background.
    radiative_part = offset+1:offset+2count
    secondary_tortoise_rhs!(
        view(derivative, radiative_part), view(state, radiative_part), tortoise_radius,
        primary_values, spin, alpha, radiative_frequency, radiative_modes, couplings,
        polar_values, azimuthal_values, grid,
    )
    include_source || return derivative

    geometry = tortoise_geometry(tortoise_radius, spin)
    secondary_values = Dict(
        mode => state[index] for (index, mode) in enumerate(secondary_modes)
    )
    source = radiative_source(
        primary_values, secondary_values, geometry.radius, spin, alpha,
        radiative_modes, polar_values, azimuthal_values, grid,
    )

    # Convert from the r-coordinate equation to dQ/dr_*, as for the secondary.
    source_scale = geometry.delta / geometry.radial_factor^2
    for (index, mode) in enumerate(radiative_modes)
        derivative[offset+count+index] -= source_scale * source[mode]
    end

    return derivative
end

"""
    radiative_outer_boundary_state(position, spin, alpha, frequency, modes, couplings,
                          amplitudes; far_position=5000.0, far_step=0.25, outer_order=16)

Return the m = 3 field's state `[R₃; Q₃]` at tortoise `position` (the primary's
outer boundary) for a purely outgoing wave with leading amplitudes `amplitudes`
at large r, in `modes` order.

ω₃ lies just above mu, so the outgoing series converges only far out: at
r ≈ 140 its terms do not decrease. Start it at tortoise `far_position`
instead and integrate the Kerr equation inward to `position` with step
`far_step`. Beyond the primary's domain the source has decayed like
exp(-(κ₁ + 2κ₂) r) and is omitted, as is the primary's averaged background
(see the comment in the body).
"""
function radiative_outer_boundary_state(position, spin, alpha, frequency, modes, couplings,
                               amplitudes; far_position=5000.0, far_step=0.25, outer_order=16)
    all(mode -> last(mode) == 3, modes) ||
        throw(ArgumentError("radiative modes must have m = 3"))
    real(frequency) > alpha ||
        throw(DomainError(frequency, "the m = 3 field must radiate (frequency > alpha)"))

    far_state = boundary_state(
        :outer, far_position, spin, alpha, frequency, modes, couplings, amplitudes;
        fundamental_m=3, outer_order,
    )
    # The primary's averaged background is left out of this leg. Beyond the
    # primary's domain it is the tail of the primary's radiating harmonics,
    # <1 - cos ϕ₁> ≲ 1e-7, falling like 1/r^2. Over the leg it mostly shifts
    # the phase of each outer solution (by up to 1.5e-5 at A1 = 13), which only
    # relabels the outer amplitudes the solve determines anyway. The part that
    # can change the solution, outside the span of the outer solutions, is
    # ≤ 1e-7 (tests/adiabatic/check_radiative_far_background.jl). For the same
    # reason the flux to infinity is taken from the current at `position`, not
    # from the amplitudes.
    rhs! = (derivative, state, tortoise_radius) -> linear_tortoise_rhs!(
        derivative, state, tortoise_radius, spin, alpha, frequency, modes, couplings;
        fundamental_m=3,
    )
    return integrate_radial(rhs!, far_state, far_position, position; step=far_step)
end

"""
    radiative_matching_residual(horizon_amplitudes, outer_amplitudes, frequency,
                                secondary_solution, primary_profiles, primary_modes,
                                spin, alpha, secondary_modes, radiative_modes,
                                couplings, polar_values, azimuthal_values, grid;
                                include_source=true, horizon_order=3, outer_order=7,
                                far_position=5000.0, far_step=0.25, far_order=16,
                                save_profile=false)

Return the mismatch `[R_horizon - R_outer; Q_horizon - Q_outer]` of the m = 3
field at the primary's matching point.

Each side starts the secondary from `secondary_solution` (from
[`solve_secondary_mode`](@ref)). The m = 3 field starts from the ingoing
horizon series at the inner end, and from
[`radiative_outer_boundary_state`](@ref) at the outer end, with the given
amplitudes. Both are then integrated together with
[`radiative_tortoise_rhs!`](@ref) on the primary's domain at twice its step.
`horizon_order` and `outer_order` are the secondary's (and the horizon
series') orders; `far_position`, `far_step`, and `far_order` set the m = 3
field's far outer leg. Only the m = 3 part enters the residual. With
`save_profile=true`, also return both integrated branches.
"""
function radiative_matching_residual(
    horizon_amplitudes, outer_amplitudes, frequency, secondary_solution,
    primary_profiles, primary_modes, spin, alpha, secondary_modes, radiative_modes,
    couplings, polar_values, azimuthal_values, grid;
    include_source=true, horizon_order=3, outer_order=7,
    far_position=5000.0, far_step=0.25, far_order=16, save_profile=false,
)
    primary_from_horizon = primary_profiles.from_horizon
    primary_from_outer = primary_profiles.from_outer
    inner = first(primary_from_horizon.positions)
    match = last(primary_from_horizon.positions)
    outer = first(primary_from_outer.positions)

    secondary_m = last(first(secondary_modes))
    secondary_frequency = secondary_solution.frequency
    secondary_state(side, position, amplitudes) = boundary_state(
        side, position, spin, alpha, secondary_frequency, secondary_modes, couplings,
        amplitudes; fundamental_m=secondary_m, horizon_order, outer_order,
    )
    horizon_state = vcat(
        secondary_state(:horizon, inner, secondary_solution.horizon_amplitudes),
        boundary_state(
            :horizon, inner, spin, alpha, frequency, radiative_modes, couplings,
            horizon_amplitudes; fundamental_m=3, horizon_order,
        ),
    )
    outer_state = vcat(
        secondary_state(:outer, outer, secondary_solution.outer_amplitudes),
        radiative_outer_boundary_state(
            outer, spin, alpha, frequency, radiative_modes, couplings, outer_amplitudes;
            far_position, far_step, outer_order=far_order,
        ),
    )

    # Right-hand side that reads the primary from one saved branch.
    rhs_reading(primary_branch) = (derivative, state, position) ->
        radiative_tortoise_rhs!(
            derivative, state, position,
            primary_values_at(primary_branch, position, primary_modes),
            spin, alpha, secondary_frequency, frequency, secondary_modes,
            radiative_modes, couplings, polar_values, azimuthal_values, grid;
            include_source,
        )

    # Twice the primary spacing puts each RK4 half-step on a primary sample,
    # as in secondary_matching_residual.
    primary_spacing = abs(primary_from_horizon.positions[2] - primary_from_horizon.positions[1])
    step = 2 * primary_spacing * (1 + 1e-9)

    from_horizon = integrate_radial(
        rhs_reading(primary_from_horizon), horizon_state, inner, match; step, save_profile,
    )
    from_outer = integrate_radial(
        rhs_reading(primary_from_outer), outer_state, outer, match; step, save_profile,
    )

    # The m = 3 field follows the secondary's 2N entries.
    radiative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2length(radiative_modes)
    if save_profile
        residual = from_horizon.states[radiative_part, end] - from_outer.states[radiative_part, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon[radiative_part] - from_outer[radiative_part]
end

"""
    solve_radiative_mode(frequency, secondary_solution, primary_profiles,
                         primary_modes, spin, alpha, secondary_modes,
                         radiative_modes, couplings, polar_values,
                         azimuthal_values, grid; series_options...)

Solve for the m = 3 field at `frequency` (2 Re ω₂ - Re ω₁) driven by the
primary and the secondary, and return its amplitudes, profiles, and fluxes.

The residual of [`radiative_matching_residual`](@ref) is affine in the 2N
boundary amplitudes, so solve it by superposition as
[`solve_dissipative_mode`](@ref) does. Fluxes are in units `Fa² M` per unit
(2,2) outer amplitude of the secondary to the fourth power; positive means
energy lost. `infinity_flux` is the outward current at the primary's outer
boundary, which equals the flux to infinity because the far leg has no source
and a real operator. `horizon_flux` is negative when the field is superradiant
(ω₃ < 3Ω_H): the black hole then supplies energy to it.
"""
function solve_radiative_mode(frequency, secondary_solution, primary_profiles,
                              primary_modes, spin, alpha, secondary_modes,
                              radiative_modes, couplings, polar_values,
                              azimuthal_values, grid; series_options...)
    count = length(radiative_modes)
    residual_at(amplitudes; include_source) = radiative_matching_residual(
        amplitudes[1:count], amplitudes[count+1:end], frequency, secondary_solution,
        primary_profiles, primary_modes, spin, alpha, secondary_modes, radiative_modes,
        couplings, polar_values, azimuthal_values, grid; include_source, series_options...,
    )

    source_residual = residual_at(zeros(ComplexF64, 2count); include_source=true)
    matching_matrix = hcat([
        residual_at([i == j ? 1.0 + 0im : 0im for j in 1:2count]; include_source=false)
        for i in 1:2count
    ]...)
    amplitudes = -(matching_matrix \ source_residual)

    # Condition number with columns scaled to unit norm, as in solve_dissipative_mode.
    column_norms = [norm(column) for column in eachcol(matching_matrix)]
    condition_number = cond(matching_matrix ./ transpose(column_norms))

    profiles = radiative_matching_residual(
        amplitudes[1:count], amplitudes[count+1:end], frequency, secondary_solution,
        primary_profiles, primary_modes, spin, alpha, secondary_modes, radiative_modes,
        couplings, polar_values, azimuthal_values, grid; save_profile=true, series_options...,
    )

    # Loss into the horizon is minus the outward current at the inner endpoint.
    radiative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2count
    horizon_flux = -outward_secondary_flux(
        profiles.from_horizon.states[radiative_part, 1],
        first(profiles.from_horizon.positions), frequency, spin, count,
    )
    infinity_flux = outward_secondary_flux(
        profiles.from_outer.states[radiative_part, 1],
        first(profiles.from_outer.positions), frequency, spin, count,
    )

    return (;
        horizon_amplitudes=amplitudes[1:count], outer_amplitudes=amplitudes[count+1:end],
        residual_norm=norm(profiles.residual), source_residual_norm=norm(source_residual),
        condition_number, profiles, horizon_flux, infinity_flux,
    )
end
