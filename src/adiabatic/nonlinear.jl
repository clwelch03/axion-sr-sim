"""Evaluate phi-sin(phi) accurately even when phi is tiny (dimensionless field)."""
function sine_remainder(field::Real)
    if abs(field) < 0.1
        squared = field^2
        # Through phi^11. Direct subtraction loses the cubic term near zero.
        return field * squared * (1/6 + squared * (-1/120 + squared *
            (1/5040 + squared * (-1/362880 + squared / 39916800))))
    end
    return field - sin(field)
end

"""Reconstruct phi=2Re(sum R_lm Y_lm exp(i*m*varphi)) at one radius and t=0.

radial_values is a dictionary keyed by (l,m), containing only positive m.
The conjugate partners are included by 2Re; m=0 needs a separate convention.
Return a real array with polar rows and azimuthal columns. R is the paper's
unrescaled radial profile: there is no 1/r or 1/sqrt(2pi) here.
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

"""Project an angular field onto the requested (l,m) modes."""
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

"""Return the nonlinear radial source N_lm as a dictionary for the requested modes.

N_lm = mu^2 int dx dvarphi/(2pi) Y_lm exp(-i*m*varphi)
                         * (r^2+a^2*x^2) * (phi-sin(phi)).
The mass term is already in radial.jl: (Box-mu^2)phi+mu^2(phi-sin(phi))=0.
Thus N enters the r-coordinate equation with a PLUS sign; Q' receives
-source_scale*N. No Delta or tortoise-coordinate factor is included here.

The printed V'(phi) term in PDF Eq. (21) appears inconsistent with its linear
mass terms. This implementation follows the field equation, PDF Eq. (5).
modes may include harmonics absent from radial_values, to measure generation.
The full sine is retained; angular resolution must be checked at large amplitude.
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

"""Evaluate the coupled nonlinear radial RHS for state=[R; Q]."""
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