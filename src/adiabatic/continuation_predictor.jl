"""
    predict_nonlinear_mode(continuation, amplitude)

Predict the frequency and physical boundary amplitudes at a larger `amplitude`.

Supply accepted entries with `amplitude` and `solution`, as saved by the
continuation driver. Extrapolate linearly in `A0` from the last two entries,
using their actual amplitude spacing. With only one entry, retain its
frequency and scale its boundary amplitudes by the amplitude ratio.

Return a named tuple with `frequency`, `horizon_amplitudes`, and
`outer_amplitudes`. Set the first outer amplitude exactly to the requested
positive real `A0`. Require finite, increasing positive amplitudes and
compatible boundary-vector lengths. Leave the accepted solutions unchanged.
The prediction is an initial guess; it must be corrected and checked by
[`solve_nonlinear_mode`](@ref) before being accepted.

See also [`solve_nonlinear_mode`](@ref).
"""
function predict_nonlinear_mode(continuation, amplitude)
    isempty(continuation) && throw(ArgumentError("require an accepted solution"))
    previous = last(continuation)
    isfinite(amplitude) && isfinite(previous.amplitude) &&
        0 < previous.amplitude < amplitude ||
        throw(ArgumentError("require a finite target A0 greater than the last positive A0"))
    previous_solution = previous.solution
    if length(continuation) == 1
        ratio = amplitude / previous.amplitude
        frequency = previous_solution.frequency
        horizon_amplitudes = ratio .* previous_solution.horizon_amplitudes
        outer_amplitudes = ratio .* previous_solution.outer_amplitudes
    else
        earlier = continuation[end-1]
        isfinite(earlier.amplitude) && 0 < earlier.amplitude < previous.amplitude ||
            throw(ArgumentError("require increasing positive continuation amplitudes"))
        earlier_solution = earlier.solution
        length(earlier_solution.horizon_amplitudes) == length(previous_solution.horizon_amplitudes) &&
            length(earlier_solution.outer_amplitudes) == length(previous_solution.outer_amplitudes) ||
            throw(DimensionMismatch("require matching boundary-vector lengths"))
        extrapolation = (amplitude - previous.amplitude) /
            (previous.amplitude - earlier.amplitude)
        frequency = previous_solution.frequency +
            extrapolation * (previous_solution.frequency - earlier_solution.frequency)
        horizon_amplitudes = previous_solution.horizon_amplitudes .+
            extrapolation .* (previous_solution.horizon_amplitudes .- earlier_solution.horizon_amplitudes)
        outer_amplitudes = previous_solution.outer_amplitudes .+
            extrapolation .* (previous_solution.outer_amplitudes .- earlier_solution.outer_amplitudes)
    end
    isempty(outer_amplitudes) && throw(ArgumentError("require a fundamental outer amplitude"))
    outer_amplitudes[1] = amplitude
    return (; frequency, horizon_amplitudes, outer_amplitudes)
end

"""
    predict_nonlinear_mode_cubic(continuation, amplitude; history_spacing=0.05)

Predict frequency and physical boundary amplitudes by cubic extrapolation.

Select four accepted points near equally spaced offsets behind the latest
amplitude, using their actual unequal spacing in the interpolation weights.
Fall back to [`predict_nonlinear_mode`](@ref) when history is insufficient
or the selected nodes are too close. Fix the fundamental outer amplitude
exactly to `amplitude`. Leave accepted solutions unchanged.

# Arguments
- `history_spacing=0.05`: positive finite spacing for selecting history in A0.

The prediction is an initial guess requiring correction and residual checks.
"""
function predict_nonlinear_mode_cubic(continuation, amplitude; history_spacing=0.05)
    fallback = predict_nonlinear_mode(continuation, amplitude)
    isfinite(history_spacing) && history_spacing > 0 ||
        throw(ArgumentError("require positive finite history spacing"))
    nodes = [entry.amplitude for entry in continuation]
    all(value -> isfinite(value) && value > 0, nodes) &&
        all(diff(nodes) .> 0) ||
        throw(ArgumentError("require increasing positive finite history"))
    indices = [findlast(value -> value <= last(nodes) - offset * history_spacing,
                       nodes) for offset in 3:-1:0]
    any(isnothing, indices) && return fallback
    selected = continuation[Int.(indices)]
    nodes = [entry.amplitude for entry in selected]
    all(diff(nodes) .>= history_spacing / 2) || return fallback
    for entry in selected
        length(entry.solution.horizon_amplitudes) == length(fallback.horizon_amplitudes) &&
            length(entry.solution.outer_amplitudes) == length(fallback.outer_amplitudes) ||
            throw(DimensionMismatch("require matching boundary-vector lengths"))
    end
    weights = [prod((amplitude - nodes[other]) / (nodes[index] - nodes[other])
                   for other in eachindex(nodes) if other != index)
               for index in eachindex(nodes)]
    frequency = sum(weight * entry.solution.frequency
                    for (weight, entry) in zip(weights, selected))
    horizon_amplitudes = sum(weight .* entry.solution.horizon_amplitudes
                            for (weight, entry) in zip(weights, selected))
    outer_amplitudes = sum(weight .* entry.solution.outer_amplitudes
                          for (weight, entry) in zip(weights, selected))
    outer_amplitudes[1] = amplitude
    isfinite(frequency) && all(isfinite, horizon_amplitudes) &&
        all(isfinite, outer_amplitudes) || return fallback
    return (; frequency, horizon_amplitudes, outer_amplitudes)
end
