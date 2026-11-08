"""
Kerr background geometry in Boyer–Lindquist coordinates with `M = 1`.

Everything downstream is expressed in units of the black-hole mass, so the
only free geometric parameter is the dimensionless spin `a = J/M² ∈ [0, 1)`.

**Radial coordinate convention.** Near the horizon `r − r₊` is as small as
`1e-11` on the tortoise domains we use, so forming it from `r` loses most of
its digits. All radial functions therefore take the *horizon offset*
`s = r − r₊` rather than `r`; use `radius(background, s)` when `r` is needed.
"""

"""
    KerrBackground(a)

Kerr black hole of unit mass and dimensionless spin `a`. Stores the derived
horizon radii `r₊`, `r₋`, their difference `horizon_gap = r₊ − r₋`, and the
horizon angular velocity `Ω_H = a/(r₊² + a²)`.
"""
struct KerrBackground{T<:Real}
    a::T
    rplus::T
    rminus::T
    horizon_gap::T
    horizon_angular_velocity::T
    function KerrBackground(a::T) where {T<:Real}
        0 ≤ a < 1 || throw(DomainError(a, "spin must satisfy 0 ≤ a < 1"))
        root = sqrt(1 - a^2)
        rplus, rminus = 1 + root, 1 - root
        new{T}(a, rplus, rminus, rplus - rminus, a / (rplus^2 + a^2))
    end
end

"""Boyer–Lindquist radius `r = r₊ + s`."""
radius(background::KerrBackground, s) = background.rplus + s

"""Horizon function `Δ = (r − r₊)(r − r₋) = s(s + r₊ − r₋)`, exact near the horizon."""
delta(background::KerrBackground, s) = s * (s + background.horizon_gap)

"""Radial derivative `Δ' = 2r − 2`."""
delta_prime(background::KerrBackground, s) = 2radius(background, s) - 2

"""`r² + a² = 2r₊ + 2r₊ s + s²` (using `r₊² + a² = 2r₊`)."""
r2a2(background::KerrBackground, s) = 2background.rplus + 2background.rplus * s + s^2

"""`ρ² = r² + a² x²` with `x = cos θ`."""
rho2(background::KerrBackground, s, x) = radius(background, s)^2 + background.a^2 * x^2

"""Jacobian `dr/dr* = Δ/(r² + a²)`."""
dr_drstar(background::KerrBackground, s) = delta(background, s) / r2a2(background, s)

"""
    tortoise(background, s)

Tortoise coordinate `r*` with `dr*/dr = (r² + a²)/Δ`, fixing the integration
constant so that `r* = r + (2r₊/(r₊−r₋)) ln s − (2r₋/(r₊−r₋)) ln(s + r₊ − r₋)`.
Near the horizon this reduces to `r* ≈ (2r₊/(r₊−r₋)) ln s + const`.
"""
function tortoise(background::KerrBackground, s)
    gap = background.horizon_gap
    return radius(background, s) + (2background.rplus / gap) * log(s) -
           (2background.rminus / gap) * log(s + gap)
end

"""
    offset_from_tortoise(background, rstar; tol=1e-13)

Invert `tortoise`, returning `s = r − r₊`. Newton iteration runs in
`log_offset = ln s`, where the map is smooth and monotone all the way down the
throat, with a bisection safeguard inside a bracket.
"""
function offset_from_tortoise(background::KerrBackground, rstar; tol=1e-13, maxiter=200)
    mismatch(log_offset) = tortoise(background, exp(log_offset)) - rstar
    # d r*/d(ln s) = s · dr*/dr
    mismatch_slope(log_offset) = (offset = exp(log_offset); offset / dr_drstar(background, offset))

    # The horizon asymptote ln s ≈ (r* − r₊)(r₊−r₋)/(2r₊) and the far-field
    # asymptote ln s ≈ ln r* are each accurate in their own regime; widen from
    # both until the sign change is enclosed.
    horizon_estimate = (rstar - background.rplus) * background.horizon_gap / (2background.rplus)
    far_field_estimate = log(max(rstar, 1.0))
    bracket_low = min(horizon_estimate, far_field_estimate) - 1
    bracket_high = max(horizon_estimate, far_field_estimate) + 1
    while mismatch(bracket_low) > 0
        bracket_low -= 1
    end
    while mismatch(bracket_high) < 0
        bracket_high += 1
    end

    log_offset = rstar > 5 ? far_field_estimate : horizon_estimate
    for _ in 1:maxiter
        current_mismatch = mismatch(log_offset)
        converged = abs(current_mismatch) < tol * max(1, abs(rstar)) || bracket_high - bracket_low < 1e-15
        converged && return exp(log_offset)
        current_mismatch > 0 ? (bracket_high = log_offset) : (bracket_low = log_offset)
        newton_step = log_offset - current_mismatch / mismatch_slope(log_offset)
        inside_bracket = bracket_low < newton_step < bracket_high
        log_offset = inside_bracket ? newton_step : (bracket_low + bracket_high) / 2
    end
    error("offset_from_tortoise did not converge at r* = $rstar")
end
