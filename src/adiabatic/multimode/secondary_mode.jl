"""
    secondary_background_source(primary_values, secondary_values, radius, spin, alpha,
                                secondary_modes, polar_values, azimuthal_values, grid;
                                reading=:averaged)

Return the background term of Omiya et al. 2211.01949 Eq. (61) for each
secondary mode at one radius, as a dictionary keyed by mode.

With `reading=:averaged` (the default), replace `1 - cos ϕ₁` by its average over
one period of the primary, [`averaged_mass_deficit`](@ref), and project
`alpha^2 * (r^2 + a^2 cos^2θ) * <1 - cos ϕ₁> * ϕ₂` onto `secondary_modes`.
This keeps only the part of the product at the secondary's own (ω₂, m₂),
as the paper's (ω, m) classification (Sec. II B) requires.

With `reading=:fixed_time`, reconstruct ϕ₁ and ϕ₂ at t = 0 and project their
product directly. That version also picks up a term at (4ω₁ - ω₂, m₂)
proportional to conj(R), which makes ω₂ depend on the secondary's
normalization phase; it is kept only to reproduce that check.

The term enters the secondary's radial equation in the slot of the nonlinear
source `N`, including the factor mu^2 missing from the printed Eq. (61).
"""
function secondary_background_source(primary_values, secondary_values, radius, spin, alpha,
                                     secondary_modes, polar_values, azimuthal_values, grid;
                                     reading=:averaged)
    # rho^2 = r^2 + a^2 cos^2θ at each polar node.
    rho_squared = radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2

    if reading === :averaged
        secondary_m = last(first(secondary_modes))
        all(mode -> last(mode) == secondary_m, secondary_modes) ||
            throw(ArgumentError("all secondary modes must share one azimuthal number"))

        mass_deficit_average = averaged_mass_deficit(
            primary_values, polar_values, azimuthal_values, grid,
        )
        # The secondary's polar profile sum_l R_l Y_l(x) at each polar node.
        secondary_profile = sum(
            value .* polar_values[mode] for (mode, value) in secondary_values
        )
        weighted_term = alpha^2 .* grid.polar_weights .* rho_squared .*
                        mass_deficit_average .* secondary_profile

        return Dict(mode => sum(weighted_term .* polar_values[mode])
                    for mode in secondary_modes)

    elseif reading === :fixed_time
        primary_field = reconstruct_field(primary_values, polar_values, azimuthal_values, grid)
        secondary_field = reconstruct_field(secondary_values, polar_values, azimuthal_values, grid)

        # 2 sin²(ϕ₁/2) equals 1 - cos ϕ₁ without cancellation when ϕ₁ is small.
        background_term =
            alpha^2 .* rho_squared .* 2 .* sin.(primary_field ./ 2) .^ 2 .* secondary_field

        return project_field(background_term, secondary_modes, polar_values, azimuthal_values, grid)
    end

    throw(ArgumentError("reading must be :averaged or :fixed_time"))
end


"""
    secondary_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                            spin, alpha, frequency, secondary_modes, couplings,
                            polar_values, azimuthal_values, grid; reading=:averaged)

Write the radial derivative of the secondary's Eq. (61) into `derivative` and return it.

Use the state layout of [`linear_tortoise_rhs!`](@ref): `state = [R; Q]` in
`secondary_modes` order, all modes sharing one azimuthal number `m2`, with
`frequency = omega2`. Supply the primary's `R_lm` at this `tortoise_radius` as
`primary_values`. Add the background term from
[`secondary_background_source`](@ref) to the slope derivatives after tortoise
scaling, as [`nonlinear_tortoise_rhs!`](@ref) does for the primary.
"""
function secondary_tortoise_rhs!(derivative, state, tortoise_radius, primary_values,
                                 spin, alpha, frequency, secondary_modes, couplings,
                                 polar_values, azimuthal_values, grid; reading=:averaged)
    # With fundamental_m = m2, every secondary mode has n = 1 and frequency omega2.
    secondary_m = last(first(secondary_modes))
    linear_tortoise_rhs!(
        derivative, state, tortoise_radius, spin, alpha,
        frequency, secondary_modes, couplings; fundamental_m=secondary_m,
    )

    count = length(secondary_modes)
    geometry = tortoise_geometry(tortoise_radius, spin)

    # Only R enters the background term; Q contains radial derivatives.
    secondary_values = Dict(
        mode => state[index] for (index, mode) in enumerate(secondary_modes)
    )
    source = secondary_background_source(
        primary_values, secondary_values, geometry.radius, spin, alpha,
        secondary_modes, polar_values, azimuthal_values, grid; reading,
    )

    # Convert from the r-coordinate equation to dQ/dr_*, as for the primary.
    source_scale = geometry.delta / geometry.radial_factor^2
    for (index, mode) in enumerate(secondary_modes)
        derivative[count + index] -= source_scale * source[mode]
    end

    return derivative
end

"""
    primary_values_at(branch, position, primary_modes)

Return the primary's `R_lm` at `position` from a saved profile branch, keyed by mode.

`branch` holds `positions` and `states` as returned by
[`integrate_radial`](@ref) with `save_profile=true`. Throw if `position` is
not one of the saved samples, so a misaligned secondary step cannot silently
read the wrong radius.
"""
function primary_values_at(branch, position, primary_modes)
    primary_spacing = branch.positions[2] - branch.positions[1]
    sample = round(Int, (position - branch.positions[1]) / primary_spacing) + 1

    (1 <= sample <= length(branch.positions) &&
        abs(branch.positions[sample] - position) <= 1e-9) ||
        throw(ArgumentError("position $position is not a saved primary sample"))

    return Dict(
        mode => branch.states[k, sample] for (k, mode) in enumerate(primary_modes)
    )
end

"""
    secondary_matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                                primary_profiles, primary_modes, spin, alpha,
                                secondary_modes, couplings, polar_values,
                                azimuthal_values, grid;
                                horizon_order=3, outer_order=7, save_profile=false, reading=:averaged)

Return the mismatch `[R_horizon - R_outer; Q_horizon - Q_outer]` of the
secondary's Eq. (61) at the primary's matching point.

`primary_profiles` is the output of [`nonlinear_matching_residual`](@ref) with
`save_profile=true` for one accepted primary solution. The secondary uses the
primary's domain: its inner, matching, and outer positions. Each secondary branch
reads the primary from the primary branch covering the same interval, using
twice the primary's spacing so that every RK4 stage falls on a saved sample.
Boundary states follow Eqs. (58)-(59) with `(omega1, m1)` replaced by
`(omega2, m2)`.
"""
function secondary_matching_residual(
    frequency, horizon_amplitudes, outer_amplitudes, primary_profiles, primary_modes,
    spin, alpha, secondary_modes, couplings, polar_values, azimuthal_values, grid;
    horizon_order=3, outer_order=7, save_profile=false, reading=:averaged
)
    primary_from_horizon = primary_profiles.from_horizon
    primary_from_outer = primary_profiles.from_outer

    # The primary's domain: horizon branch runs inner -> match, outer branch outer -> match.
    inner = first(primary_from_horizon.positions)
    match = last(primary_from_horizon.positions)
    outer = first(primary_from_outer.positions)
    isapprox(last(primary_from_outer.positions), match; atol=1e-12) ||
        throw(ArgumentError("primary branches must meet at one matching point"))

    secondary_m = last(first(secondary_modes))
    horizon_state = boundary_state(
        :horizon, inner, spin, alpha, frequency, secondary_modes, couplings,
        horizon_amplitudes; fundamental_m=secondary_m, horizon_order,
    )
    outer_state = boundary_state(
        :outer, outer, spin, alpha, frequency, secondary_modes, couplings,
        outer_amplitudes; fundamental_m=secondary_m, outer_order,
    )

    # Right-hand side that reads the primary from one saved branch.
    rhs_reading(primary_branch) = (derivative, state, position) ->
        secondary_tortoise_rhs!(
            derivative, state, position,
            primary_values_at(primary_branch, position, primary_modes),
            spin, alpha, frequency, secondary_modes, couplings,
            polar_values, azimuthal_values, grid,
        )

    # Twice the primary spacing puts each RK4 half-step on a primary sample.
    # The tiny enlargement keeps integrate_radial's ceil() from adding a step
    # through floating-point rounding; primary_values_at throws if alignment fails.
    primary_spacing = abs(primary_from_horizon.positions[2] - primary_from_horizon.positions[1])
    step = 2 * primary_spacing * (1 + 1e-9)

    from_horizon = integrate_radial(
        rhs_reading(primary_from_horizon), horizon_state, inner, match; step, save_profile,
    )
    from_outer = integrate_radial(
        rhs_reading(primary_from_outer), outer_state, outer, match; step, save_profile,
    )

    if save_profile
        residual = from_horizon.states[:, end] - from_outer.states[:, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon - from_outer
end

"""
    solve_secondary_mode(frequency_guess, primary_profiles, primary_modes, spin, alpha,
                         secondary_modes, couplings, polar_values, azimuthal_values, grid;
                         horizon_guess, outer_tail_guess, tolerance=1e-8,
                         max_iterations=15, difference_step=1e-6,
                         jacobian_backend=:finite, series_options...)

Solve the secondary's Eq. (61) on a fixed primary for `omega2` and the boundary amplitudes.

Fix the first outer amplitude (the `(2, 2)` mode; put it first in
`secondary_modes`) to one. This sets the secondary's scale and, because the
fixed-time background term involves `conj(R)`, also its phase relative to the
primary. Solve with a real Jacobian ([`real_newton`](@ref)) for the frequency,
all horizon amplitudes, and the remaining outer amplitudes.

The default tolerance is looser than the primary's: with unit outer amplitude
the matching residual has a round-off floor near 1e-10.
`series_options...` forwards `horizon_order` and `outer_order` to
[`secondary_matching_residual`](@ref).
"""
function solve_secondary_mode(
    frequency_guess, primary_profiles, primary_modes, spin, alpha,
    secondary_modes, couplings, polar_values, azimuthal_values, grid;
    horizon_guess, outer_tail_guess, tolerance=1e-8, max_iterations=15,
    difference_step=1e-6, jacobian_backend=:finite, reading=:averaged,
    series_options...,
)
    count = length(secondary_modes)
    length(horizon_guess) == count && length(outer_tail_guess) == count - 1 ||
        throw(DimensionMismatch("require N horizon and N-1 outer guesses"))

    # Unknowns: frequency, N horizon amplitudes, N-1 free outer amplitudes.
    initial = vcat(
        complex(frequency_guess), complex.(horizon_guess), complex.(outer_tail_guess),
    )
    residual = state -> secondary_matching_residual(
        state[1], state[2:count+1], vcat(one(state[1]), state[count+2:end]),
        primary_profiles, primary_modes, spin, alpha, secondary_modes, couplings,
        polar_values, azimuthal_values, grid; reading, series_options...,
    )

    # Keep the secondary on the positive-frequency bound branch.
    admissible = state -> all(isfinite, state) && 0 < real(state[1]) < alpha

    solution = real_newton(
        residual, initial;
        tolerance, max_iterations, difference_step, jacobian_backend, admissible,
    )

    return (;
        frequency=solution.state[1],
        horizon_amplitudes=solution.state[2:count+1],
        outer_amplitudes=vcat(one(solution.state[1]), solution.state[count+2:end]),
        residual=solution.residual,
        residual_norm=norm(solution.residual),
        iterations=solution.iterations,
    )
end
