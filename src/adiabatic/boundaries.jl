# Evaluate the asymptotic series from boundary_series.jl at finite endpoints.
# Units and time convention: M=1, exp(-i*omega*t). frequency is omega_n
# for the actual harmonic, whereas linear_tortoise_rhs! accepts omega_0.
# The outputs are vectors in l_values order, ready to form state=[R; Q].
# They are boundary guesses for a future matching solve, not a solved cloud.

"""Choose the decaying (:bound) or positive-frequency outgoing (:outgoing) sheet.

An explicit choice keeps a bound mode decaying even for Im(omega)<0.
Exactly zero wave number needs a different expansion and is rejected.
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

"""Evaluate a vector polynomial and its derivative with respect to its argument.

For y(z)=v_0+v_1*z+...+v_N*z^N, Horner's rule repeatedly forms z*y+v_j,
starting with v_N. Differentiating that update gives z*y'+y. z is s near
the horizon and 1/r at the outer boundary; the caller applies the chain rule.
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

"""Ingoing boundary values (R,Q=dR/dr_*) for one same-frequency, same-m block.

amplitudes are coefficients of s^p, s=r-r_plus (a local normalization).
Linear Frobenius corrections are included through the requested order.
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

"""Outer boundary values including inverse-radius corrections (PDF Eq. 23).

amplitudes are the coefficients multiplying exp(i*k*r)*r^p at infinity;
the fundamental bound coefficient is A0 in PDF Eq. 20, with M=1.
order=7 retains terms through 1/r^7. Finite endpoint/order convergence must
still be checked in a profile solve. These are linear asymptotic boundaries.
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

"""Assemble boundary state=[R; Q] for all modes, preserving their ordering.

frequency is omega_0; amplitudes follows the same ordering as modes.
side is :horizon or :outer. Boundary series remain linear asymptotic
approximations, whose accuracy must be checked by moving the endpoints.
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

    state = zeros(ComplexF64, 2count)

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