"""
    sine_remainder!(output::AbstractArray{Float64}, field::AbstractArray{Float64})

Evaluate `field - sin(field)` elementwise, write the values into `output`, and return it.

Require matching array axes. Use a degree-eleven Taylor polynomial for
`abs(field[index]) < 0.1` to avoid cancellation, and direct subtraction otherwise.
Use the same threshold and polynomial as `nonlinear_source.py`, without fast math.
The input is unchanged unless it aliases `output`; using the same array for
both arguments is supported, but overlapping views with different index mappings
are not supported.

See also [`sine_remainder`](@ref).
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

"""
    sine_remainder(field::AbstractArray{Float64})

Return an array containing the elementwise remainder `field - sin(field)`.

Allocate the result with `similar(field)` and leave `field` unchanged.
Use the stable evaluation and threshold of [`sine_remainder!`](@ref).
"""
function sine_remainder(field::AbstractArray{Float64})
    return sine_remainder!(similar(field), field)
end
