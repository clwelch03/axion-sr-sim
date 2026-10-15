# Coordinate utilities only: these functions do not evolve the scalar field.
# Units: M=1, so spin=a/M is also the numerical value of a.
# r is the Boyer-Lindquist radius; r_* is the tortoise integration coordinate.
# Both describe the same exterior point, with dr/dr_* = Delta/(r^2+a^2).

"""
    kerr_horizons(spin)

Return the Kerr horizon radii and their separation in units `M = 1`.

Require finite `spin = a/M` with `abs(spin) < 1`. The returned named tuple
contains `r_plus`, `r_minus`, and `gap = r_plus - r_minus`. Evaluate the inner
radius without subtractive cancellation at small spin. Extremal Kerr is excluded.

# Examples
```jldoctest
julia> kerr_horizons(0.0)
(r_plus = 2.0, r_minus = 0.0, gap = 2.0)
```
"""
function kerr_horizons(spin)
    isfinite(spin) && abs(spin) < 1 || throw(DomainError(spin, "require |spin| < 1"))
    half_gap = sqrt((1 - spin) * (1 + spin))
    r_plus = 1 + half_gap
    # r_+ * r_- = a^2. This avoids subtracting two almost equal numbers in
    # r_- = 1-sqrt(1-a^2) when the spin is small.
    r_minus = spin^2 / r_plus
    return (; r_plus, r_minus, gap=2half_gap)
end

# Using z=log(r-r_plus) preserves the horizon displacement even when adding
# it to r_plus rounds away. The logarithms use M=1; this fixes the additive
# constant in r_* = r + (2r_+/gap)log(r-r_+) - (2r_-/gap)log(r-r_-).
"""
    tortoise_from_log_distance(log_distance, horizons)

Evaluate the tortoise coordinate from `log_distance = log(r - r_plus)`.

Supply `horizons` from [`kerr_horizons`](@ref). In units `M = 1`, fix the additive
constant by `r_* = r + (2r_plus/gap)*log(r-r_plus) - (2r_minus/gap)*log(r-r_minus)`.
The logarithmic displacement preserves points whose radius rounds to `r_plus`.

See also [`tortoise_coordinate`](@ref), [`tortoise_geometry`](@ref).
"""
function tortoise_from_log_distance(log_distance, horizons)
    distance = exp(log_distance)
    # s=distance, g=gap. Compute log(s+g) by factoring out the larger term.
    # log1p(x) evaluates log(1+x) accurately even when x is very small.
    log_outer_distance = log_distance > log(horizons.gap) ?
        log_distance + log1p(horizons.gap / distance) :
        log(horizons.gap) + log1p(distance / horizons.gap)
    # Since r_+-r_-=g, the logarithmic part can be regrouped as
    # 2log(s+g) + (2r_+/g)[log(s)-log(s+g)]. This is the same r_*(r)
    # relation as above, with one explicit choice of additive constant.
    return horizons.r_plus + distance + 2log_outer_distance +
           (2horizons.r_plus / horizons.gap) * (log_distance - log_outer_distance)
end

"""
    tortoise_coordinate(horizon_distance, spin)

Map a positive horizon displacement `horizon_distance = r - r_plus` to `r_*`.

Require a finite, positive displacement and finite `abs(spin) < 1`.
Use units `M = 1` and the additive constant defined by
[`tortoise_from_log_distance`](@ref).

See also [`tortoise_geometry`](@ref).
"""
function tortoise_coordinate(horizon_distance, spin)
    isfinite(horizon_distance) && horizon_distance > 0 ||
        throw(DomainError(horizon_distance, "require positive finite horizon distance"))
    return tortoise_from_log_distance(log(horizon_distance), kerr_horizons(spin))
end

"""
    tortoise_geometry(tortoise_radius, spin)

Find the exterior Kerr geometry at the finite tortoise coordinate `tortoise_radius`.

Solve for `log(r - r_plus)` with finite `abs(spin) < 1` in units `M = 1`.
Return a named tuple containing `radius`, `delta = Δ`, `radial_factor = r^2+a^2`,
`dr_drstar = Δ/(r^2+a^2)`, `horizon_distance`, `log_horizon_distance`,
`r_plus`, `r_minus`, and `gap`.

Near the horizon, `radius` can round to `r_plus` while the displacement remains
positive. Use the displacement fields to preserve that information.
Throw an error if inversion fails to converge or the displacement underflows.

See also [`tortoise_coordinate`](@ref).
"""
function tortoise_geometry(tortoise_radius, spin)
    isfinite(tortoise_radius) || throw(DomainError(tortoise_radius, "require finite r_*"))
    horizons = kerr_horizons(spin)
    # Seed Newton's method with the near-horizon logarithmic approximation
    # on the inner side, or the large-radius approximation r_* ~ r outside.
    # These are only starting guesses; the full mapping determines the result.
    log_distance = tortoise_radius < horizons.r_plus ?
        (tortoise_radius - horizons.r_plus +
         2*horizons.r_minus / horizons.gap * log(horizons.gap)) *
        horizons.gap / (2*horizons.r_plus) :
        log(max(tortoise_radius - horizons.r_plus, one(tortoise_radius)))
    tolerance = 16 * eps(typeof(float(tortoise_radius + spin))) * max(1, abs(tortoise_radius))

    # Stop on the error in r_*, scaled to its floating-point magnitude.
    for iteration in 1:100
        distance = exp(log_distance)
        radius = horizons.r_plus + distance
        radial_factor = radius^2 + spin^2
        residual = tortoise_from_log_distance(log_distance, horizons) - tortoise_radius
        if abs(residual) <= tolerance
            distance > 0 || throw(DomainError(tortoise_radius, "horizon distance underflow"))
            # Delta=(r-r_+)(r-r_-)=s(s+g). The polynomial r^2-2r+a^2
            # loses accuracy here because its large terms nearly cancel.
            delta = distance * (distance + horizons.gap)
            return (; radius, delta, radial_factor, dr_drstar=delta / radial_factor,
                      horizon_distance=distance, log_horizon_distance=log_distance,
                      horizons...)
        end
        # d r_*/d log(r-r_+) = (r^2+a^2)/(r-r_-).
        correction = residual * (distance + horizons.gap) / radial_factor
        # Limit the change in z during an iteration to avoid a large jump
        # from an imperfect initial guess. This is a numerical safeguard,
        # not an approximation to the coordinate relation.
        log_distance -= clamp(correction, -2, 2)
    end
    error("tortoise inversion did not converge")
end

"""
    kerr_metric(tortoise_radius, spin, cosine_theta)

Returns the value of the Kerr metric at a given point. 
"""
function kerr_metric(tortoise_radius, spin, cosine_theta)
    -1 < cosine_theta < 1 ||
        throw(DomainError(cosine_theta, "require -1 < cos(theta) < 1"))

    geometry = tortoise_geometry(tortoise_radius, spin)
    radius, delta = geometry.radius, geometry.delta
    sine_squared = 1 - cosine_theta^2
    rho_squared = radius^2 + spin^2 * cosine_theta^2
    azimuth_factor = geometry.radial_factor^2 - spin^2 * delta * sine_squared

    metric = zeros(typeof(float(rho_squared)), 4, 4)
    metric[1, 1] = -(delta - spin^2 * sine_squared) / rho_squared
    metric[1, 4] = metric[4, 1] = -2spin * radius * sine_squared / rho_squared
    metric[2, 2] = rho_squared / delta
    metric[3, 3] = rho_squared
    metric[4, 4] = azimuth_factor * sine_squared / rho_squared

    inverse_metric = zero(metric)
    inverse_metric[1, 1] = -azimuth_factor / (rho_squared * delta)
    inverse_metric[1, 4] = inverse_metric[4, 1] =
        -2spin * radius / (rho_squared * delta)
    inverse_metric[2, 2] = delta / rho_squared
    inverse_metric[3, 3] = 1 / rho_squared
    inverse_metric[4, 4] = (delta - spin^2 * sine_squared) /
        (rho_squared * delta * sine_squared)

    return (; metric, inverse_metric)
end