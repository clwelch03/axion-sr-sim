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
