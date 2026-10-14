"""
    component_frequency_and_m(q, cloud_frequencies, cloud_m)

Calculate the frequency and summed azimuthal integer for a Fourier component.

Return the named tuple `(frequency, m)`, with
`frequency = sum(q .* cloud_frequencies)` and `m = sum(q .* cloud_m)`.

# Example
```julia
numbers = component_frequency_and_m((2, -1), [0.405, 0.415], [1, 2])
# numbers.frequency ≈ 0.395; numbers.m == 0
```
"""
function component_frequency_and_m(q, cloud_frequencies, cloud_m)
    length(q) == length(cloud_frequencies) == length(cloud_m) ||
        throw(DimensionMismatch("one q entry per independent cloud"))

    frequency = sum(q .* cloud_frequencies)
    m = sum(q .* cloud_m)
    return (; frequency, m)
end

"""
    fourier_component(q, l_values, cloud_m)

Define a retained Fourier component and validate its angular indices.

Supply an integer vector or tuple `q` listing the independent clouds, and
their respective azimuthal numbers `cloud_m`. The Fourier component's azimuthal
number is then sum(q .* cloud_m). Also supply integer `l_values`.

Return `(q=Tuple(q), l_values=collect(l_values))`: an immutable label tuple
and a copied vector of retained angular indices. No frequency or azimuthal
number is stored. Use [`component_frequency_and_m`](@ref) to calculate them
for given cloud frequencies.

# Example
```julia
component = fourier_component((2, -1), [0, 2, 4], [1, 2])
# Retain l = 0, 2, 4 for the driven component with m = 0.
```
"""
function fourier_component(q, l_values, cloud_m)
    length(q) == length(cloud_m) ||
        throw(DimensionMismatch("need one q entry per independent cloud"))
    !isempty(q) && all(value -> value isa Integer, q) &&
        all(value -> value isa Integer, cloud_m) ||
        throw(ArgumentError("q and cloud_m must contain integers and be nonempty"))

    m = sum(q .* cloud_m)
    !isempty(l_values) && all(l -> l isa Integer && l >= abs(m), l_values) ||
        throw(ArgumentError("require nonempty integer l_values with l >= |m|"))
    allunique(l_values) ||
        throw(ArgumentError("l_values must be distinct"))

    return (; q=Tuple(q), l_values=collect(l_values))
end

"""
    cloud_label(n::Integer, l::Integer, m::Integer)

Label an independent cloud by its linear-limit quantum numbers.

Here `n` is the principal quantum number. Require `n > l >= abs(m)`,
and return the named tuple `(n, l, m)`.

This records the intended state; it does not calculate its profile,
frequency, or amplitude normalization.
"""
function cloud_label(n::Integer, l::Integer, m::Integer)
    n > l >= abs(m) ||
        throw(ArgumentError("require n > l >= |m|"))
    return (; n, l, m)
end
