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

Compute the coefficients of the radial equation in Boyer–Lindquist radius.

Use units `M = 1`, with `spin = a/M`, `alpha = M*μ`, and `mode_frequency`
equal to the actual harmonic frequency, rather than the fundamental frequency.
For the unrescaled radial profile, the equation is
`(Δ*R′)′ + B*R + coupling_scale*(C*R) + N = 0`, with primes denoting `d/dr`.

Return a named tuple with `delta = Δ`, `delta_prime = dΔ/dr`,
`local_coefficient = B`, and `coupling_scale = spin^2*(mode_frequency^2-alpha^2)`.
The local coefficient excludes all angular coupling terms, including the
diagonal of `C`. This is the radial operator in PDF Eq. (21).

See also [`tortoise_operator_coefficients`](@ref).
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

Compute the polar matrix elements of `cos(θ)^2` by quadrature.

Return a dictionary keyed by `(l, other_l, m)` for `0 ≤ m ≤ l, other_l ≤ l_max`.
Each entry integrates `Y_lm(x) * x^2 * Y_other_l,m(x)` over `x = cos(θ)` in
`[-1, 1]`. Include diagonal entries; entries outside `other_l = l, l ± 2`
vanish up to quadrature roundoff.

Use tables from [`build_polar_values`](@ref) and weights from [`angular_grid`](@ref).
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
"""
    radial_second_derivative(radial_value, radial_slope, coupled_value,
                             nonlinear_source, coefficients)

Return `d²R/dr²` by rearranging the radial equation for one mode.

Supply `radial_slope = dR/dr`, `coupled_value = (C*R)[l,m]` including the
diagonal term, and the source `N` before tortoise scaling as `nonlinear_source`.
Obtain `coefficients` from [`radial_operator_coefficients`](@ref).
This expression divides by `Δ`; use [`tortoise_operator_coefficients`](@ref)
to avoid that division near the horizon.
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
    tortoise_operator_coefficients(geometry, spin, alpha, mode_frequency, l, m)

Compute radial coefficients that remain regular near the horizon in tortoise coordinates.

Use `geometry` from [`tortoise_geometry`](@ref), units `M = 1`, and the actual
harmonic frequency `mode_frequency`. Return a named tuple defining
`Q′ = -slope_coefficient*Q - local_coefficient*R - coupling_scale*(C*R) - source_scale*N`,
where `Q = dR/dr_*` and the prime means `d/dr_*`. The source `N` belongs to the
original radial equation. Combine factors analytically to avoid division by `Δ`.
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
