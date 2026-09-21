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
    reconstruct_field(radial_values, polar_values, azimuthal_values, grid)

Reconstruct the real angular field at one radius and time zero.

Supply `radial_values` keyed by `(l, m)` with `0 < m ≤ l`, and angular tables
from [`build_polar_values`](@ref) and [`build_azimuthal_values`](@ref).
Return a real matrix with polar rows and azimuthal columns on `grid`.

Evaluate `ϕ = 2*real(sum(R_lm * Y_lm * exp(im*m*φ)))`; the factor of two
includes the conjugate partners. The radial profile is unrescaled, with no
`1/r` or `1/sqrt(2π)` factor. Zero and negative `m` are not supported.

See also [`project_field`](@ref).
"""
function reconstruct_field(radial_values, polar_values, azimuthal_values, grid)
    # Rows are polar nodes; columns are azimuthal nodes.
    field = zeros(length(grid.cosine_theta_nodes), length(grid.azimuth_nodes))

    for ((l, m), amplitude) in radial_values
        0 < m <= l || throw(ArgumentError("require positive m <= l"))

        # transpose preserves exp(+imϕ); ' would conjugate it.
        harmonic = polar_values[(l, m)] * transpose(azimuthal_values[m])

        # Include the complex-conjugate partner to make the field real.
        field .+= 2 .* real.(amplitude .* harmonic)
    end

    return field
end

"""
    project_field(field, modes, polar_values, azimuthal_values, grid)

Project an angular field onto the requested `(l, m)` modes.

Supply `field` as a matrix with polar rows and azimuthal columns on `grid`.
Return a dictionary keyed by each requested mode, integrating against the
conjugate harmonic with measure `d(cos(θ)) * dφ/(2π)`.

See also [`reconstruct_field`](@ref), [`nonlinear_source`](@ref).
"""
function project_field(field, modes, polar_values, azimuthal_values, grid)
    # Integrate over dx dϕ/(2π). Polar weights broadcast down the rows.
    weighted_field = field .* grid.polar_weights .* grid.azimuth_weight

    return Dict(
        (l, m) => sum(
            weighted_field .* polar_values[(l, m)] .*
            transpose(conj.(azimuthal_values[m]))
        )
        for (l, m) in modes
    )
end

"""
    nonlinear_source(radial_values, radius, spin, alpha, polar_values,
                     azimuthal_values, grid; modes=keys(radial_values))

Return the nonlinear radial source as a dictionary keyed by the requested modes.

Reconstruct the field from `radial_values` and project
`alpha^2 * (radius^2 + spin^2*cos(θ)^2) * (ϕ - sin(ϕ))` onto `modes`, using
the measure and normalization of [`project_field`](@ref). Require `0 < m ≤ l`;
`modes` may include harmonics absent from `radial_values` to measure generation.
Use units `M = 1`, so `alpha` is the scalar mass in these units.

The source `N` enters the radial equation with a plus sign and contributes
`-source_scale*N` to `dQ/dr_*`. No tortoise factor is included here.
This follows the field equation in PDF Eq. (5); the printed potential term in
Eq. (21) appears inconsistent with its linear mass terms. Retain the full sine,
with a small-field evaluation from [`sine_remainder`](@ref), and check angular
resolution at large amplitude.
"""
function nonlinear_source(radial_values, radius, spin, alpha,
                          polar_values, azimuthal_values, grid;
                          modes=keys(radial_values))
    all(mode -> 0 < mode[2] <= mode[1], modes) ||
        throw(ArgumentError("require positive m <= l for projected modes"))
    field = reconstruct_field(radial_values, polar_values, azimuthal_values, grid)
    # Polar vectors broadcast along rows. azimuth_weight integrates dvarphi/(2pi).
    rho_squared = radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2
    weighted_source = alpha^2 .* rho_squared .* sine_remainder.(field) .*
                      grid.polar_weights .* grid.azimuth_weight
    # Projection uses exp(-i*m*varphi), unlike reconstruction above.
    return Dict((l, m) => sum(weighted_source .* polar_values[(l, m)] .*
                             transpose(conj.(azimuthal_values[m])))
                for (l, m) in modes)
end

"""
    nonlinear_tortoise_rhs!(derivative, state, tortoise_radius, spin, alpha,
                            frequency, modes, couplings, polar_values,
                            azimuthal_values, grid; fundamental_m=1)

Write the coupled nonlinear radial derivative into `derivative` and return it.

Use the state layout, harmonic frequencies, units, and nonaliasing requirements
of [`linear_tortoise_rhs!`](@ref). Leave `state` unchanged and add the source
from [`nonlinear_source`](@ref) to the slope derivatives after tortoise scaling.
Supply angular tables evaluated on `grid` for every retained mode.
"""
function nonlinear_tortoise_rhs!(
    derivative, state, tortoise_radius, spin, alpha, frequency,
    modes, couplings, polar_values, azimuthal_values, grid;
    fundamental_m=1,
    )
    # Start with the linear equations, including the harmonic frequencies.
    linear_tortoise_rhs!(
        derivative, state, tortoise_radius, spin, alpha,
        frequency, modes, couplings; fundamental_m,
    )

    count = length(modes)
    geometry = tortoise_geometry(tortoise_radius, spin)

    # Only R enters the local potential; Q contains radial derivatives.
    radial_values = Dict(
        mode => state[index] for (index, mode) in enumerate(modes)
    )

    source = nonlinear_source(
        radial_values, geometry.radius, spin, alpha,
        polar_values, azimuthal_values, grid; modes,
    )

    # Convert the source from the r-coordinate equation to dQ/dr_*.
    source_scale = geometry.delta / geometry.radial_factor^2
    for (index, mode) in enumerate(modes)
        derivative[count + index] -= source_scale * source[mode]
    end

    return derivative
end