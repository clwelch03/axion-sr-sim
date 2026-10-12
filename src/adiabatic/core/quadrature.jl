"""
    integrate_uniform(values, spacing)

Integrate uniformly spaced samples with composite Simpson's rule. With an odd
number of intervals, the last three use Simpson's 3/8 rule. Both rules are
exact for cubics.
"""
function integrate_uniform(values, spacing)
    interval_count = length(values) - 1
    interval_count >= 2 || throw(ArgumentError("require at least two intervals"))

    simpson_end = iseven(interval_count) ? interval_count : interval_count - 3
    total = zero(float(eltype(values)))
    if simpson_end >= 2
        part = values[1:simpson_end+1]
        total += spacing / 3 *
            (part[1] + part[end] + 4sum(part[2:2:end-1]) + 2sum(part[3:2:end-2]))
    end
    if !iseven(interval_count)
        part = values[end-3:end]
        total += 3spacing / 8 * (part[1] + 3part[2] + 3part[3] + part[4])
    end
    return total
end


"""
    cumulative_trapezoid(times, values)

Return the cumulative integral of sampled `values` over `times`, starting
at zero at the first sample. Both inputs must be nonempty vectors of equal
length, with finite entries and strictly increasing times. Return a vector
of Float64 values.
"""
function cumulative_trapezoid(times, values)
    length(times) == length(values) ||
        throw(DimensionMismatch("require one value per time"))
    isempty(times) && throw(ArgumentError("require at least one sample"))
    all(isfinite, times) && all(isfinite, values) ||
        throw(ArgumentError("times and values must be finite"))
    all(diff(times) .> 0) ||
        throw(ArgumentError("times must increase strictly"))

    integral = zeros(Float64, length(times))
    for index in 2:length(times)
        spacing = times[index] - times[index - 1]
        integral[index] = integral[index - 1] +
            spacing * (values[index - 1] + values[index]) / 2
    end
    return integral
end
