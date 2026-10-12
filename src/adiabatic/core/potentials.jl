"""
    sine_remainder(field::Real)

Evaluate `field - sin(field)` accurately for a dimensionless scalar field.

For `abs(field) < 0.1`, use the Taylor polynomial through degree eleven to
avoid cancellation; otherwise subtract directly.

# Examples
```jldoctest
julia> sine_remainder(0.0)
0.0
```
"""
function sine_remainder(field::Real)
    if abs(field) < 0.1
        squared = field^2
        # Through phi^11. Direct subtraction loses the cubic term near zero.
        return field * squared * (1/6 + squared * (-1/120 + squared *
            (1/5040 + squared * (-1/362880 + squared / 39916800))))
    end
    return field - sin(field)
end

"""
    potential_value(field::Real, potential::Symbol=:cosine)

Evaluate the dimensionless scalar potential `U(field)`.

Select `:cosine`, `:quartic`, or `:sextic`. Exclude the overall factor
`alpha^2`, which the energy calculation supplies separately. Evaluate
the cosine potential without small-field cancellation.

See also [`potential_remainder`](@ref).
"""
function potential_value(field::Real, potential::Symbol=:cosine)
    if potential === :cosine
        return 2sin(field / 2)^2
    elseif potential === :quartic
        return field^2 / 2 - field^4 / 24
    elseif potential === :sextic
        return field^2 / 2 - field^4 / 24 + field^6 / 720
    end

    throw(ArgumentError("potential must be :cosine, :quartic, or :sextic"))
end

"""
    potential_remainder(field::Real, potential::Symbol=:cosine)

Evaluate `field - U'(field)` for the selected scalar potential.

Return the nonlinear remainder after separating the linear mass term.
Exclude the factor `alpha^2` and the geometric factor supplied by
[`nonlinear_source`](@ref). Use [`sine_remainder`](@ref) for accurate
evaluation of the cosine case near zero.

See also [`potential_value`](@ref).
"""
function potential_remainder(field::Real, potential::Symbol=:cosine)
    if potential === :cosine
        return sine_remainder(field)
    elseif potential === :quartic
        return field^3 / 6
    elseif potential === :sextic
        return field^3 / 6 - field^5 / 120
    end

    throw(ArgumentError("potential must be :cosine, :quartic, or :sextic"))
end
