"""
    interpolate_branch_energy(amplitude, left, right)

Interpolate energy and its amplitude derivative between two branch samples.

Supply endpoints with fields `amplitude`, `energy`, and `energy_derivative`.
Require increasing endpoint amplitudes and an evaluation amplitude inside
their closed interval. Supply finite real energies and derivatives with
consistent normalization.

Construct a cubic Hermite interpolant matching both endpoint energies
and derivatives. Return a named tuple with `energy` and `energy_derivative`.
The latter is the derivative of the same interpolating polynomial.

Check interpolation accuracy separately from the endpoint derivative accuracy.
"""
function interpolate_branch_energy(amplitude, left, right)
    spacing = right.amplitude - left.amplitude
    spacing > 0 ||
        throw(ArgumentError("require increasing endpoint amplitudes"))

    left.amplitude <= amplitude <= right.amplitude ||
        throw(DomainError(amplitude, "amplitude lies outside the interval"))

    coordinate = (amplitude - left.amplitude) / spacing
    energy_change = right.energy - left.energy

    linear_coefficient = spacing * left.energy_derivative
    quadratic_coefficient = 3energy_change -
        spacing * (2left.energy_derivative + right.energy_derivative)
    cubic_coefficient = -2energy_change +
        spacing * (left.energy_derivative + right.energy_derivative)

    energy = left.energy + coordinate * (
        linear_coefficient + coordinate * (
            quadratic_coefficient + coordinate * cubic_coefficient
        )
    )

    energy_derivative = (
        linear_coefficient + coordinate * (
            2quadratic_coefficient + 3coordinate * cubic_coefficient
        )
    ) / spacing

    return (; energy, energy_derivative)
end

"""
    branch_interval(amplitude, amplitudes)

Find the saved branch interval containing `amplitude`.

Supply at least two finite, strictly increasing amplitudes. Require the
evaluation amplitude to lie within the saved range; do not extrapolate.

Return a named tuple containing `left_index` and `right_index`.
At an interior sample, select the interval to its right. At the final
sample, select the last interval.
"""
function branch_interval(amplitude, amplitudes)
    length(amplitudes) >= 2 ||
        throw(ArgumentError("require at least two branch samples"))

    all(isfinite, amplitudes) && all(diff(amplitudes) .> 0) ||
        throw(ArgumentError("require finite, strictly increasing amplitudes"))

    first(amplitudes) <= amplitude <= last(amplitudes) ||
        throw(DomainError(amplitude, "amplitude lies outside the saved branch"))

    left_index = min(
        searchsortedlast(amplitudes, amplitude),
        length(amplitudes) - 1,
    )

    return (; left_index, right_index=left_index + 1)
end

"""
    amplitude_rate_on_interval(amplitude, left, right)

Evaluate dA0/d(t/M) = -total_flux / energy_derivative.

Endpoints must contain amplitude, energy, energy_derivative, horizon_flux,
outer_flux, horizon_flux_derivative, and outer_flux_derivative. All derivatives
are with respect to amplitude. Energy and flux share a consistent normalization;
positive flux denotes energy lost by the cloud.

Use cubic Hermite interpolation for energy and each signed flux.
Reject an exactly zero energy derivative.
"""
function amplitude_rate_on_interval(amplitude, left, right)
    interpolated = interpolate_branch_energy(amplitude, left, right)
    derivative = interpolated.energy_derivative
    iszero(derivative) &&
        throw(DomainError(amplitude, "energy derivative is zero"))

    function interpolate_flux(value_field, derivative_field)
        left_flux = (; amplitude=left.amplitude,
            energy=getproperty(left, value_field),
            energy_derivative=getproperty(left, derivative_field))
        right_flux = (; amplitude=right.amplitude,
            energy=getproperty(right, value_field),
            energy_derivative=getproperty(right, derivative_field))
        return interpolate_branch_energy(amplitude, left_flux, right_flux).energy
    end

    horizon_flux = interpolate_flux(:horizon_flux, :horizon_flux_derivative)
    outer_flux = interpolate_flux(:outer_flux, :outer_flux_derivative)
    return -(horizon_flux + outer_flux) / derivative
end

"""
    amplitude_rate(amplitude, branch)

Evaluate the amplitude growth rate on a tabulated stationary branch.

Supply `branch.amplitudes` as the sorted amplitude grid and
`branch.samples` as the corresponding endpoint records accepted by
[`amplitude_rate_on_interval`](@ref). Each record's amplitude must match
its grid entry.

Select the containing interval and return `dA0/d(t/M)`.
Reject evaluation outside the saved branch.

See also [`branch_interval`](@ref), [`amplitude_rate_on_interval`](@ref).
"""
function amplitude_rate(amplitude, branch)
    length(branch.amplitudes) == length(branch.samples) ||
        throw(DimensionMismatch("require one sample per amplitude"))

    interval = branch_interval(amplitude, branch.amplitudes)

    return amplitude_rate_on_interval(
        amplitude,
        branch.samples[interval.left_index],
        branch.samples[interval.right_index],
    )
end

"""
    interpolate_branch_state(amplitude, left, right)

Construct a Newton guess by cubic Hermite interpolation of the normalized
stationary state. Endpoints contain `amplitude`, `solution`, and a packed-real
`state_derivative` of (frequency, horizon/A0, outer_tail/A0). Return physical
frequency and boundary amplitudes, fixing the first outer amplitude to A0.
The guess must be corrected by the stationary solver before use as a solution.
"""
function interpolate_branch_state(amplitude, left, right)
    0 < left.amplitude < right.amplitude ||
        throw(ArgumentError("require positive increasing endpoint amplitudes"))
    left.amplitude <= amplitude <= right.amplitude ||
        throw(DomainError(amplitude, "amplitude lies outside the interval"))
    mode_count = length(left.solution.horizon_amplitudes)
    for entry in (left, right)
        length(entry.solution.horizon_amplitudes) == mode_count &&
            length(entry.solution.outer_amplitudes) == mode_count &&
            length(entry.state_derivative) == 4mode_count ||
            throw(DimensionMismatch("endpoint states and tangents must use the same modes"))
    end
    normalized_state(entry) = vcat(entry.solution.frequency,
        entry.solution.horizon_amplitudes ./ entry.amplitude,
        entry.solution.outer_amplitudes[2:end] ./ entry.amplitude)
    spacing = right.amplitude - left.amplitude
    fraction = (amplitude - left.amplitude) / spacing
    state = (2fraction^3 - 3fraction^2 + 1) .* normalized_state(left) .+
        (fraction^3 - 2fraction^2 + fraction) * spacing .* unpack_complex(left.state_derivative) .+
        (-2fraction^3 + 3fraction^2) .* normalized_state(right) .+
        (fraction^3 - fraction^2) * spacing .* unpack_complex(right.state_derivative)
    return (; frequency=state[1],
        horizon_amplitudes=amplitude .* state[2:mode_count+1],
        outer_amplitudes=amplitude .* vcat(one(state[1]), state[mode_count+2:end]))
end
