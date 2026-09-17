"""
    sine_remainder!(output, field)

Evaluate Phi-sin(Phi) for real Float64 arrays, writing into output.
The ! denotes a function that modifies an argument (here, output).
Same threshold and polynomial as nonlinear_source.py; no fast-math.
"""
function sine_remainder!(output::AbstractArray{Float64}, field::AbstractArray{Float64})
    axes(output) == axes(field) || throw(DimensionMismatch("output and field must match"))
    for index in eachindex(output, field)
        value = field[index]
        if abs(value) < 0.1
            squared = value * value
            polynomial = squared / 39916800 - 1 / 362880
            polynomial = polynomial * squared + 1 / 5040
            polynomial = polynomial * squared - 1 / 120
            polynomial = polynomial * squared + 1 / 6
            output[index] = (value * squared) * polynomial
        else
            output[index] = value - sin(value)
        end
    end
    return output
end

"""Allocating wrapper, matching the Python function's input/output behavior."""
function sine_remainder(field::AbstractArray{Float64})
    return sine_remainder!(similar(field), field)
end
