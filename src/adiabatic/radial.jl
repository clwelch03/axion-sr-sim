# Radial equation conventions (M=1):
#   (Delta*R_r)_r + B*R + a^2*(omega^2-mu^2)*(C*R) + N = 0.
# R is the paper's R_tilde: there is no extra 1/r field rescaling here.
# N denotes the projected nonlinear remainder, with a plus sign in this
# equation. It is an input to the scalar helper; the assembled RHS is linear.
# mode_frequency is the actual harmonic frequency n*omega_0, not omega_0.
# Complex frequencies are squared algebraically; abs2(omega) would change
# the eigenvalue equation. Equation numbers refer to the paper's PDF.

"""
    radial_operator_coefficients(radius, spin, alpha, mode_frequency, l, m)

Calculate the factors that multiply the radial profiles and their
derivatives in the differential equation for R_{nl}.

This is PDF Eq. (21) with the diagonal and off-diagonal angular corrections
collected into C. local_coefficient is B and excludes every C contribution.
Primes in this function's interpretation are derivatives with respect to r.
"""
function radial_operator_coefficients(radius, spin, alpha, mode_frequency, l, m)
    delta = radius^2 - 2*radius + spin^2
    delta_prime = 2*radius - 2
    # K=(r^2+a^2)*omega-a*m: the frequency/rotation combination in Eq. (21).
    frequency_combination = (radius^2 + spin^2) * mode_frequency - spin * m

    local_coefficient = frequency_combination^2 / delta - alpha^2 * radius^2 +
                        2*spin * m * mode_frequency -
                        spin^2 * mode_frequency^2 - l * (l + 1)

    coupling_scale = spin^2 * (mode_frequency^2 - alpha^2)
    return (; delta, delta_prime, local_coefficient, coupling_scale)
end

"""
    build_cosine_squared_couplings(l_max, polar_values, grid)

Build a grid of angular couplings
C_{l,l'}^m = int_{-1}^1 Y_{lm}(x) x^2 Y_{l',m}(x) dx.
"""
function build_cosine_squared_couplings(l_max, polar_values, grid)
    # x=cos(theta), so the polar integration measure is dx. The azimuthal
    # projection already enforces equal m; no azimuthal weights enter here.
    weighted_cosine_squared = grid.polar_weights .* grid.cosine_theta_nodes.^2

    # Key order is (output l, input l', fixed m). The diagonal is included.
    # Analytically only l'=l,l+-2 survive. Other stored entries should be
    # zero up to quadrature roundoff, which our tests can check.
    return Dict(
        (l, other_l, m) => sum(
            weighted_cosine_squared .*
            polar_values[(l, m)] .* polar_values[(other_l, m)]
        )
        for m in 0:l_max for l in m:l_max for other_l in m:l_max
    )
end
"""Return d^2R/dr^2 for one mode by rearranging the r-coordinate equation.

radial_slope=dR/dr. coupled_value=sum(C[l,l',m]*R[l',m]) includes the
diagonal term. nonlinear_source is N before any tortoise-coordinate scaling.
This helper contains divisions by Delta; the tortoise helper below combines
those factors analytically before evaluating them near the horizon.
"""
function radial_second_derivative(
    radial_value, radial_slope, coupled_value, nonlinear_source, coefficients
)
    return -(
        coefficients.delta_prime * radial_slope +
        coefficients.local_coefficient * radial_value +
        coefficients.coupling_scale * coupled_value +
        nonlinear_source
    ) / coefficients.delta
end

"""
Coefficients of Q' = -slope_coefficient*Q - local_coefficient*R
                   - coupling_scale*(C*R) - source_scale*N, with Q=dR/dr_*.
Primes here mean d/dr_*. N is the source in the original r-coordinate
equation. No division by Delta is evaluated, including near the horizon.
"""
function tortoise_operator_coefficients(geometry, spin, alpha, mode_frequency, l, m)
    radius, delta, radial_factor = geometry.radius, geometry.delta, geometry.radial_factor
    # Write A=r^2+a^2 and Q=dR/dr_*. Then Delta*R_r=A*Q and
    # (Delta*R_r)_r = 2r*Q + (A^2/Delta)*dQ/dr_*.
    # Multiplying the whole radial equation by Delta/A^2 isolates dQ/dr_*.
    source_scale = delta / radial_factor^2
    frequency_combination = radial_factor * mode_frequency - spin * m
    regular_term = -alpha^2 * radius^2 + 2spin * m * mode_frequency -
                   spin^2 * mode_frequency^2 - l * (l + 1)
    # B=K^2/Delta+regular_term. Simplify (Delta/A^2)*B to the following
    # expression before evaluation: it has no artificial 1/Delta singularity.
    local_coefficient = (frequency_combination / radial_factor)^2 + source_scale * regular_term
    return (; slope_coefficient=2radius * source_scale, local_coefficient,
              coupling_scale=source_scale * spin^2 * (mode_frequency^2 - alpha^2),
              source_scale)
end

"""
Linear first-order system in r_* for state=[R; Q], Q=dR/dr_*.
Modes are (l,m), with positive m=n*fundamental_m; frequency is omega_0.
Couplings use (l,other_l,m) keys. derivative and state must not alias.
This computes an RHS only; it neither integrates nor solves for frequency.

For N retained modes, state[1:N] stores R and state[N+1:2N] stores Q,
in the same order as modes. A complex state stores phases as well as amplitudes.
The trailing ! means derivative is modified in place; state is read only.
"""
function linear_tortoise_rhs!(derivative, state, tortoise_radius, spin, alpha,
                              frequency, modes, couplings; fundamental_m=1)
    count = length(modes)
    length(state) == length(derivative) == 2count || throw(DimensionMismatch("expected [R; Q]"))
    fundamental_m isa Integer && fundamental_m > 0 || throw(ArgumentError("require positive fundamental_m"))
    geometry = tortoise_geometry(tortoise_radius, spin)
    for (index, (l, m)) in enumerate(modes)
        m > 0 && m % fundamental_m == 0 || throw(ArgumentError("require m=n*fundamental_m, n>=1"))
        # All harmonics share one pattern: n=m/m0 and omega_n=n*omega0.
        # Here frequency is omega0; the coefficient helper takes omega_n.
        coefficients = tortoise_operator_coefficients(
            geometry, spin, alpha, (m ÷ fundamental_m) * frequency, l, m)
        coupled_value = zero(eltype(state))
        for (other_index, (other_l, other_m)) in enumerate(modes)
            # Linear Kerr angular mixing preserves m. Nonlinear interactions
            # between different m blocks will enter through the future source.
            other_m == m || continue
            coupled_value += couplings[(l, other_l, m)] * state[other_index]
        end
        # First half: dR/dr_*=Q. Second half: the transformed radial equation.
        # No nonlinear source is included in this linear validation system.
        derivative[index] = state[count + index]
        derivative[count + index] = -coefficients.slope_coefficient * state[count + index] -
            coefficients.local_coefficient * state[index] - coefficients.coupling_scale * coupled_value
    end
    return derivative
end
