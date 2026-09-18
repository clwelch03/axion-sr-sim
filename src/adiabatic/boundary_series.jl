using LinearAlgebra: I, Diagonal

# These functions derive LINEAR asymptotic series for a single (omega,m)
# block containing several l values. They prepare boundary data; they do not
# find an eigenfrequency or integrate the profile through the cloud.
# frequency is the actual harmonic frequency n*omega_0, in M=1 units.
# amplitudes supplies the leading coefficient for each l in l_values.
# Each coefficients[j+1] is a vector v_j in that same l order: Julia arrays
# start at 1, whereas the series powers start at 0.
#
# Multiplying the radial equation by Delta gives a polynomial-coefficient ODE:
#   Delta^2 R'' + Delta*Delta' R' + (K^2 I + Delta D) R = 0,
#   K=(r^2+a^2)*omega-a*m,  D=-mu^2*r^2*I + constant_potential.
# The primes in this file mean d/dr. We substitute a series into this ODE
# and require each successive power of the expansion variable to vanish.

# Linear angular block in Delta*R'' + Delta'*R' + (K^2/Delta + D)*R=0,
# with D=-alpha^2*r^2*I + constant_potential. All entries have the same m, omega.
function boundary_angular_block(spin, alpha, frequency, l_values, m, couplings, amplitudes)
    count = length(l_values)
    count > 0 && length(amplitudes) == count || throw(DimensionMismatch("one amplitude per l"))
    length(unique(l_values)) == count && all(l -> l >= abs(m), l_values) ||
        throw(ArgumentError("require distinct l >= |m|"))
    identity_matrix = Matrix{typeof(complex(frequency + alpha + spin))}(I, count, count)
    # The dictionary becomes a matrix for this one fixed-m block. Its row
    # and column order must both follow l_values to act on amplitude vectors.
    angular_matrix = [couplings[(l, other_l, m)] for l in l_values, other_l in l_values]
    constant_potential = (2spin * m * frequency - spin^2 * frequency^2) * identity_matrix -
        Diagonal([l * (l + 1) for l in l_values]) +
        spin^2 * (frequency^2 - alpha^2) * angular_matrix
    return identity_matrix, constant_potential
end

"""Return the ingoing exponent p and vectors v_0,...,v_order near the horizon.

The ansatz R=s^p sum(v_j*s^j), s=r-r_plus, factors out the non-polynomial
wave behavior before expanding the remaining smooth part in powers of s.
This is a Frobenius series. v_0=amplitudes fixes the local normalization;
it need not equal the horizon amplitude convention in the prefactor of
PDF Eq. (22). The ingoing exponent is the same.
"""
function horizon_series(spin, alpha, frequency, l_values, m, couplings, amplitudes; order=3)
    order isa Integer && order >= 0 || throw(ArgumentError("require nonnegative series order"))
    horizons = kerr_horizons(spin)
    identity_matrix, constant = boundary_angular_block(
        spin, alpha, frequency, l_values, m, couplings, amplitudes)
    radius, gap = horizons.r_plus, horizons.gap
    # Expand K about r_+: K=k0+k1*s+k2*s^2.
    k0, k1, k2 = (radius^2 + spin^2) * frequency - spin * m, 2radius * frequency, frequency
    # The lowest-order equation is (gap^2*p^2+k0^2)*v_0=0.
    # Choose p=-i*k0/gap: because r_*~(2r_+/gap)log(s), this is the
    # ingoing factor exp[-i*(omega-m*Omega_H)*r_*], up to a constant.
    exponent = -im * k0 / gap
    # Likewise D=d0+d1*s+d2*s^2. These are matrices across l, not scalars.
    d0 = constant - alpha^2 * radius^2 * identity_matrix
    d1 = -2alpha^2 * radius * identity_matrix
    d2 = -alpha^2 * identity_matrix
    # Coefficients of K^2 + Delta*D in ascending powers of s.
    potential = [k0^2 * identity_matrix,
        2k0 * k1 * identity_matrix + gap * d0,
        (k1^2 + 2k0 * k2) * identity_matrix + d0 + gap * d1,
        2k1 * k2 * identity_matrix + d1 + gap * d2,
        k2^2 * identity_matrix + d2]
    # Delta^2 starts at s^2; Delta*Delta' starts at s.
    second = [gap^2, 2gap, one(gap)]
    first = [gap^2, 3gap, 2one(gap)]
    # Explicitly: Delta^2=s^2*(g^2+2g*s+s^2),
    # and Delta*Delta'=s*(g^2+3g*s+2s^2), where g=gap.
    # Integer input amplitudes must still allow complex floating-point corrections.
    coefficient_type = promote_type(eltype(constant), eltype(amplitudes))
    coefficients = [coefficient_type.(amplitudes)]
    for index in 1:order
        # index=j: collect the already-known contribution to s^(p+j).
        # An offset h selects v_(j-h), whose own exponent is p+j-h.
        residual = zero(coefficients[1])
        for offset in 1:min(index, 4)
            power = exponent + index - offset
            # Differentiating s^power contributes power*(power-1) or power.
            # Derivative-polynomial corrections extend only through h=2;
            # the potential extends through h=4 and mixes l components.
            scalar = offset <= 2 ?
                second[offset + 1] * power * (power - 1) + first[offset + 1] * power : 0
            previous = coefficients[index - offset + 1]
            residual += scalar * previous + potential[offset + 1] * previous
        end
        # The remaining term multiplies the unknown v_j. Using the indicial
        # equation above, its coefficient simplifies to g^2*j*(j+2p).
        # Thus denominator*v_j + residual = 0 fixes the next vector.
        denominator = gap^2 * index * (index + 2exponent)
        iszero(denominator) && throw(DomainError(exponent, "resonant horizon series"))
        push!(coefficients, -residual / denominator)
    end
    return (; exponent, coefficients)
end

# Fetch the coefficient of r^degree, treating absent powers as zero. This
# lets the outer recurrence ask for degrees -2,-1,... without special arrays.
polynomial_coefficient(coefficients, degree) =
    0 <= degree < length(coefficients) ? coefficients[degree + 1] : 0

"""Return k, p, and vectors v_0,...,v_order for the large-radius expansion.

R=exp(i*k*r)*r^p sum(v_j/r^j), with v_0=amplitudes. The power p already
includes the geometric 1/r falloff. kind selects the branch of k through
outer_wave_number in boundaries.jl. This is the form of PDF Eq. (23),
with angular mixing retained as matrix operations on each v_j.
"""
function outer_series(spin, alpha, frequency, l_values, m, couplings, amplitudes;
                      kind=:bound, order=7)
    order isa Integer && order >= 0 || throw(ArgumentError("require nonnegative series order"))
    wave_number = outer_wave_number(frequency, alpha; kind)
    identity_matrix, constant = boundary_angular_block(
        spin, alpha, frequency, l_values, m, couplings, amplitudes)
    # The leading r^4 balance requires k^2=omega^2-mu^2. The next r^3
    # balance fixes p=-1+i*(2omega^2-mu^2)/k, the long-range correction.
    exponent = -1 + im * (2frequency^2 - alpha^2) / wave_number
    exponential_rate = im * wave_number
    k0 = spin^2 * frequency - spin * m
    # Delta^2, Delta*Delta', and K^2+Delta*D, in ascending powers of r.
    second = [spin^4, -4spin^2, 4 + 2spin^2, -4, 1]
    first = [-2spin^2, 4 + 2spin^2, -6, 2]
    potential = [k0^2 * identity_matrix + spin^2 * constant,
        -2constant, (2frequency * k0 - spin^2 * alpha^2) * identity_matrix + constant,
        2alpha^2 * identity_matrix, (frequency^2 - alpha^2) * identity_matrix]
    # Use the operator's precision as well as the supplied amplitudes' type.
    coefficient_type = promote_type(eltype(constant), eltype(amplitudes))
    coefficients = [coefficient_type.(amplitudes)]
    for index in 1:order
        # At order j, match the coefficient of exp(i*k*r)*r^(p+3-j).
        # The j-th unknown has coefficient -2*i*k*j; lower series orders
        # supply the residual that we now calculate.
        residual = zero(coefficients[1])
        for offset in 1:min(index, 5)
            # v_(j-h) has a factor r^e with e=p-(j-h). To contribute to
            # the target power, the ODE must raise its power by d=3-h.
            # d can reach -2 because taking two derivatives lowers by two.
            degree = 3 - offset
            power = exponent - (index - offset)
            # For a term exp(sigma*r)*r^e, differentiation contributes
            # sigma^2, 2sigma*e/r, e(e-1)/r^2 to R'', and sigma,e/r to R'.
            # Combine those with the polynomial coefficients of Delta^2
            # and Delta*Delta'. The potential matrix is added separately.
            scalar = exponential_rate^2 * polynomial_coefficient(second, degree) +
                exponential_rate * (2power * polynomial_coefficient(second, degree + 1) +
                                    polynomial_coefficient(first, degree)) +
                power * (power - 1) * polynomial_coefficient(second, degree + 2) +
                power * polynomial_coefficient(first, degree + 1)
            previous = coefficients[index - offset + 1]
            residual += scalar * previous
            degree >= 0 && (residual += potential[degree + 1] * previous)
        end
        # -2*i*k*j*v_j + residual = 0. Exact k=0 was rejected above;
        # the threshold requires a different asymptotic expansion.
        push!(coefficients, residual / (2exponential_rate * index))
    end
    return (; wave_number, exponent, coefficients)
end
