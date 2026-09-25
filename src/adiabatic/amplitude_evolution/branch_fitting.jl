"""
    local_energy_derivative(amplitudes, energies, center_amplitude; degree=2)

Estimate the energy derivative with respect to amplitude at `center_amplitude`.

Supply equal-length vectors of amplitudes A0 and corresponding energies
from one stationary solution branch. Select the local sample neighborhood
before calling; all supplied samples enter the fit with equal weight.
Samples need not be sorted or uniformly spaced. Supply finite real values
and at least `degree + 1` distinct amplitudes.

Fit a polynomial by least squares in the centered, scaled coordinate
`x = (A0 - center_amplitude) / amplitude_scale`, where `amplitude_scale`
is the largest absolute amplitude offset. Evaluate the derivative at
`x = 0` by dividing the fitted linear coefficient by `amplitude_scale`.

Return the scalar estimate of `dE/dA0` at `center_amplitude`, preserving
the supplied energy normalization. For energies stored as `E/Fa²`,
the result is `d(E/Fa²)/dA0`.

Check sensitivity to the sample neighborhood and polynomial degree before
using the derivative in amplitude evolution. The function does not estimate
derivative uncertainty or restrict evaluation to the sampled amplitude range.

# Arguments
- `degree=2`: positive integer polynomial degree; use `2` for a quadratic
  fit or `3` for a cubic fit.

See also [`integrate_energy_branch`](@ref).
"""
function local_energy_derivative(
    amplitudes, energies, center_amplitude; degree=2,
)
    length(amplitudes) == length(energies) ||
        throw(DimensionMismatch("require one energy per amplitude"))

    degree >= 1 ||
    throw(ArgumentError("require degree >= 1"))

    length(unique(amplitudes)) >= degree + 1 ||
        throw(ArgumentError("require at least degree + 1 distinct amplitudes"))

    offsets = amplitudes .- center_amplitude
    amplitude_scale = maximum(abs, offsets)
    scaled_offsets = offsets ./ amplitude_scale
    
    design_matrix = hcat(
        [scaled_offsets.^power for power in 0:degree]...
    )

    coefficients = design_matrix \ energies
    return coefficients[2] / amplitude_scale
end