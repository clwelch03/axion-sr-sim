# Evaluate the asymptotic series from boundary_series.jl at finite endpoints.
# Units and time convention: M=1, exp(-i*omega*t). frequency is omega_n
# for the actual harmonic, whereas linear_tortoise_rhs! accepts omega_0.
# The outputs are vectors in l_values order, ready to form state=[R; Q].
# They supply asymptotic boundary data to the matching solvers.

"""
    outer_wave_number(frequency, alpha; kind=:bound)

Choose the decaying or positive-frequency outgoing wave number.

Use units `M = 1` and time dependence `exp(-im*frequency*t)`. Require finite
`real(frequency) > 0` and finite `alpha ≥ 0`. For `kind = :bound`, require
`real(frequency) < alpha` and return `im*sqrt(complex(alpha^2-frequency^2))`;
this branch decays even for negative imaginary frequency. For `kind = :outgoing`,
require `real(frequency) ≥ alpha` and return `sqrt(complex(frequency^2-alpha^2))`.
Reject zero wave number, which requires a separate threshold expansion.

See also [`outer_series`](@ref).
"""
function outer_wave_number(frequency, alpha; kind=:bound)
    isfinite(frequency) && real(frequency) > 0 || throw(DomainError(frequency, "require Re(omega)>0"))
    isfinite(alpha) && alpha >= 0 || throw(DomainError(alpha, "require alpha>=0"))
    if kind == :bound
        real(frequency) < alpha || throw(DomainError(frequency, "bound branch requires Re(omega)<mu"))
        # Write q=sqrt(mu^2-omega^2), with Re(q)>0, and k=i*q.
        # Then exp(i*k*r)=exp(-q*r) decays. Expressing this via q keeps
        # the decaying choice when the imaginary part of omega changes sign.
        wave_number = im * sqrt(complex(alpha^2 - frequency^2))
    elseif kind == :outgoing
        real(frequency) >= alpha || throw(DomainError(frequency, "outgoing branch requires Re(omega)>=mu"))
        # For positive real omega above the mass threshold, k>0 and
        # exp(i*k*r-i*omega*t) is outgoing. Continue that root for complex
        # omega; an outgoing damped mode need not decay spatially at infinity.
        wave_number = sqrt(complex(frequency^2 - alpha^2))
    else
        throw(ArgumentError("kind must be :bound or :outgoing"))
    end
    iszero(wave_number) && throw(DomainError(frequency, "zero wave number needs a threshold expansion"))
    return wave_number
end

"""
    boundary_polynomial(coefficients, argument) -> (value, derivative)

Evaluate a vector polynomial and its derivative using Horner's rule.

Store the coefficient of power `j` in `coefficients[j+1]`, starting with the
constant term. Return the value and derivative with respect to `argument`.
For boundary series, this argument is `r-r_plus` or `1/r`; apply the chain rule
separately to obtain a radial derivative.
"""
function boundary_polynomial(coefficients, argument)
    value = copy(coefficients[end])
    derivative = zero(value)
    for index in (length(coefficients) - 1):-1:1
        # Update the derivative first: it needs the previous polynomial value.
        derivative = derivative * argument + value
        value = value * argument + coefficients[index]
    end
    return value, derivative
end

"""
    horizon_boundary_data(tortoise_radius, spin, alpha, frequency, l_values,
                          m, couplings, amplitudes; order=3)

Evaluate ingoing boundary profiles and tortoise slopes for one fixed-frequency block.

Return a named tuple with vectors `values = R` and `tortoise_slopes = dR/dr_*`
in `l_values` order. The `frequency` is the actual harmonic frequency and
`amplitudes` contains the leading coefficients of `(r-r_plus)^p` in units `M = 1`.
Include linear Frobenius corrections through the nonnegative integer `order`.

See also [`horizon_series`](@ref), [`outer_boundary_data`](@ref).
"""
function horizon_boundary_data(tortoise_radius, spin, alpha, frequency,
                               l_values, m, couplings, amplitudes; order=3)
    # construct geometry and boundary polynomial
    geometry = tortoise_geometry(tortoise_radius, spin)
    series = horizon_series(spin, alpha, frequency, l_values, m, couplings, amplitudes; order)
    polynomial, derivative = boundary_polynomial(series.coefficients, geometry.horizon_distance)

    # R=s^p*y(s). exp(p*log(s)) handles the complex power using positive s.
    carrier = exp(series.exponent * geometry.log_horizon_distance)
    values = carrier * polynomial
    
    # dR/dr_* = (Delta/A)*s^p*(p*y/s + dy/ds), A=r^2+a^2.
    # Cancel the s before evaluation using Delta=s*(s+gap). This avoids
    # constructing a huge 1/s near the horizon and multiplying it by tiny Delta.
    tortoise_slopes = carrier * (geometry.horizon_distance + geometry.gap) / geometry.radial_factor *
        (series.exponent * polynomial + geometry.horizon_distance * derivative)
    return (; values, tortoise_slopes)
end

"""
    outer_boundary_data(tortoise_radius, spin, alpha, frequency, l_values,
                        m, couplings, amplitudes; kind=:bound, order=7)

Evaluate outer boundary profiles and tortoise slopes with inverse-radius corrections.

Return a named tuple with vectors `values = R` and `tortoise_slopes = dR/dr_*`
in `l_values` order. Supply the actual harmonic `frequency` in units `M = 1`.
The `amplitudes` multiply `exp(im*k*r)*r^p` at infinity; the fundamental bound
coefficient is `A0` in PDF Eq. (20).

Select `kind` as in [`outer_wave_number`](@ref). The nonnegative integer `order`
retains terms through `1/r^order`. These are linear asymptotic boundaries;
check endpoint and series-order convergence when solving a profile.

See also [`outer_series`](@ref), [`horizon_boundary_data`](@ref).
"""
function outer_boundary_data(tortoise_radius, spin, alpha, frequency,
                             l_values, m, couplings, amplitudes; kind=:bound, order=7)
    geometry = tortoise_geometry(tortoise_radius, spin)
    series = outer_series(spin, alpha, frequency, l_values, m, couplings, amplitudes; kind, order)
    radius = geometry.radius
    polynomial, derivative = boundary_polynomial(series.coefficients, inv(radius))
    # R=E(r)*y(1/r), E=exp(i*k*r)*r^p. The polynomial derivative returned
    # above is dy/d(1/r), not dy/dr.
    carrier = exp(im * series.wave_number * radius + series.exponent * log(radius))
    values = carrier * polynomial
    # Product and chain rules: E'/E=i*k+p/r, and d(1/r)/dr=-1/r^2.
    # Finally multiply by dr/dr_* to return the slope in the integration
    # coordinate, in the same convention as the horizon boundary data.
    tortoise_slopes = geometry.dr_drstar * carrier *
        ((im * series.wave_number + series.exponent / radius) * polynomial - derivative / radius^2)
    return (; values, tortoise_slopes)
end

"""
    boundary_state(side, position, spin, alpha, frequency, modes, couplings,
                   amplitudes; fundamental_m=1, horizon_order=3, outer_order=7)

Assemble the complex boundary state `[R; Q]` in `modes` order, with `Q = dR/dr_*`.

Choose `side = :horizon` or `:outer` and a tortoise `position`. Supply one
amplitude per `(l, m)` mode. The `frequency` is fundamental; each positive
`m = n*fundamental_m` receives frequency `n*frequency`, with positive integer
`fundamental_m`. At the outer boundary, select the bound or outgoing branch
according to each harmonic's real frequency relative to `alpha`.

Use nonnegative integer `horizon_order` and `outer_order` for the respective
linear asymptotic series. Check accuracy by varying orders and endpoints.

See also [`horizon_boundary_data`](@ref), [`outer_boundary_data`](@ref).
"""
function boundary_state(
    side, position, spin, alpha, frequency, modes, couplings, amplitudes;
    fundamental_m=1, horizon_order=3, outer_order=7,
)
    side in (:horizon, :outer) ||
        throw(ArgumentError("side must be :horizon or :outer"))
    fundamental_m isa Integer && fundamental_m > 0 ||
        throw(ArgumentError("require positive fundamental_m"))

    count = length(modes)
    length(amplitudes) == count ||
        throw(DimensionMismatch("require one amplitude per mode"))

    state_type = promote_type(typeof(complex(frequency)), eltype(amplitudes))
    state = zeros(state_type, 2count)

    for m in unique([m for (l, m) in modes])
        m > 0 && m % fundamental_m == 0 ||
            throw(ArgumentError("require m=n*fundamental_m, n>=1"))

        # Gather all angular modes sharing this harmonic frequency.
        indices = findall(mode -> mode[2] == m, modes)
        l_values = [modes[index][1] for index in indices]
        mode_frequency = (m ÷ fundamental_m) * frequency

        if side == :horizon
            data = horizon_boundary_data(
                position, spin, alpha, mode_frequency,
                l_values, m, couplings, amplitudes[indices];
                order=horizon_order,
            )
        else
            kind = real(mode_frequency) < alpha ? :bound : :outgoing
            data = outer_boundary_data(
                position, spin, alpha, mode_frequency,
                l_values, m, couplings, amplitudes[indices];
                kind, order=outer_order,
            )
        end

        # Scatter each block back into the full state.
        state[indices] = data.values
        state[count .+ indices] = data.tortoise_slopes
    end

    return state
end
