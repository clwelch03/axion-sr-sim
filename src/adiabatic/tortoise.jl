# Coordinate utilities only: these functions do not evolve the scalar field.
# Units: M=1, so spin=a/M is also the numerical value of a.
# r is the Boyer-Lindquist radius; r_* is the tortoise integration coordinate.
# Both describe the same exterior point, with dr/dr_* = Delta/(r^2+a^2).

"""Return the event horizon r_plus, inner horizon r_minus, and their separation.

The logarithmic coordinate formula below applies to nonextremal Kerr (|a|<1).
Exactly extremal Kerr requires a different expression.
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

"""Map the positive horizon displacement r-r_plus to tortoise radius r_*."""
function tortoise_coordinate(horizon_distance, spin)
    isfinite(horizon_distance) && horizon_distance > 0 ||
        throw(DomainError(horizon_distance, "require positive finite horizon distance"))
    return tortoise_from_log_distance(log(horizon_distance), kerr_horizons(spin))
end

"""Find the exterior geometry at a requested tortoise coordinate.

The numerical unknown is z=log(r-r_plus), not r itself. Near the horizon,
r can round to r_plus even though the positive displacement remains
representable. Keep horizon_distance and log_horizon_distance for this reason.

The returned named tuple also contains radius=r, delta=Delta,
radial_factor=r^2+a^2, dr_drstar=dr/dr_*, and the horizon radii/gap.
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
