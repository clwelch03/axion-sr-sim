"""
    reconstruct_two_phase_field_at_radii(radial_values, cloud_m, polar_values,
                                        phase1_nodes, phase2_nodes)

Reconstruct the real field at one radius for two independent cloud phases.

Supply complex radial coefficients keyed by (l, q), with nonzero integer
2-tuples q and one representative per conjugate pair. `cloud_m` lists the two
cloud azimuthal integers.

Return a real array with axes (polar node, first phase, second phase), using
2*real(sum(R_lq * Y_l,abs(m_q) * exp(-im*q⋅psi))).
"""
function reconstruct_two_phase_field_at_radii(radial_values, cloud_m, polar_values,
                                     phase1_nodes, phase2_nodes)
    length(cloud_m) == 2 || throw(DimensionMismatch("require two cloud m values"))
    field_type = promote_type(Float64,
        (typeof(real(value)) for value in values(radial_values))...)
    field = zeros(field_type, length(first(values(polar_values))),
        length(phase1_nodes), length(phase2_nodes))

    for ((l, q), radial_value) in radial_values
        length(q) == 2 || throw(DimensionMismatch("require two entries in q"))
        m = sum(q .* cloud_m)
        polar = polar_values[(l, abs(m))]
        for j in eachindex(phase1_nodes), k in eachindex(phase2_nodes)
            phase_factor = cis(-q[1]*phase1_nodes[j] - q[2]*phase2_nodes[k])
            for i in eachindex(polar)
                field[i, j, k] += 2real(radial_value * polar[i] * phase_factor)
            end
        end
    end
    return field
end

"""
    project_two_phase_field(field, output_pairs, cloud_m, polar_values,
                            polar_weights, phase1_nodes, phase2_nodes)

Extract coefficients of requested `(l, q)` angular/phase patterns at one radius.

`field` has axes (polar node, first phase, second phase). Supply integer
2-tuples `q` and the same cloud ordering and polar table as for reconstruction.

Project onto each requested angular/phase basis function by summing field * Y_l,
abs(m_q) * exp(+im*q⋅psi) over all grid nodes, weighted by the polar quadrature
weights and divided by the number of phase pairs. Return the resulting coefficients
in a dictionary keyed by the requested `(l, q)` pairs without an extra factor of 2.
"""
function project_two_phase_field(field, output_pairs, cloud_m, polar_values,
                                 polar_weights, phase1_nodes, phase2_nodes)
    length(cloud_m) == 2 || throw(DimensionMismatch("require two cloud m values"))
    !isempty(polar_weights) && !isempty(phase1_nodes) && !isempty(phase2_nodes) ||
        throw(ArgumentError("quadrature grids must be nonempty"))
    size(field) == (length(polar_weights), length(phase1_nodes), length(phase2_nodes)) ||
        throw(DimensionMismatch("field axes must match the quadrature grids"))
    # Uniform phase quadrature is an average: divide by the number of phase pairs.
    phase_count = length(phase1_nodes) * length(phase2_nodes)

    # Each requested (l, q) becomes one dictionary entry; the begin block computes its value.
    return Dict((l, q) => begin
        length(q) == 2 || throw(DimensionMismatch("require two entries in q"))
        m = sum(q .* cloud_m)
        polar = polar_values[(l, abs(m))]
        length(polar) == length(polar_weights) ||
            throw(DimensionMismatch("polar table and weights must match"))
        # i labels cos(theta), j and k the two independent phases: field[i,j,k] is G_ijk.
        # Weights times Y integrate the polar dependence; exp(+i q⋅psi) selects q.
        sum(polar_weights[i] * polar[i] * field[i, j, k] *
            cis(q[1]*phase1_nodes[j] + q[2]*phase2_nodes[k])
            for i in eachindex(polar_weights), j in eachindex(phase1_nodes),
                k in eachindex(phase2_nodes)) / phase_count
    end for (l, q) in output_pairs)
end

"""
    nonlinear_two_phase_source(radial_values, radius, spin, alpha, output_pairs,
                               cloud_m, polar_values, grid,
                               phase1_nodes, phase2_nodes; potential=:cosine)

Compute N_(l,q) = alpha^2 * P_(l,q)[(r^2 + spin^2*x^2)*h(Phi)] at one radius.

Reconstruct the joint real field from `radial_values`, keyed by `(l, q)`.
Use the reconstruction and projector conventions for inputs and phase grids.
`radius` is Boyer-Lindquist r, not tortoise radius; x = cos(theta), `grid`
supplies `cosine_theta_nodes` and `polar_weights` matching `polar_values`.
Select `potential=:cosine`, `:quartic`, or `:sextic` for h(Phi) = Phi - U'(Phi).

Return a dictionary keyed by `output_pairs`, which may include labels absent
from the input. The radial slope equation contains
-Delta/(r^2 + spin^2)^2 times this source;
that minus sign and tortoise factor are not included here.
"""
function nonlinear_two_phase_source(radial_values, radius, spin, alpha,
                                    output_pairs, cloud_m, polar_values, grid,
                                    phase1_nodes, phase2_nodes; potential=:cosine)
    field = reconstruct_two_phase_field_at_radii(radial_values, cloud_m,
        polar_values, phase1_nodes, phase2_nodes)

    # Evaluate h on the total field so interactions between clouds are retained.
    remainder = potential_remainder.(field, Ref(potential))

    # The geometric factor varies along the polar axis, and belongs inside P.
    geometry_factor = reshape(radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2,
        :, 1, 1)

    projected = project_two_phase_field(geometry_factor .* remainder,
        output_pairs, cloud_m, polar_values, grid.polar_weights,
        phase1_nodes, phase2_nodes)

    return Dict(pair => alpha^2 * coefficient for (pair, coefficient) in projected)
end

"""
    build_polar_theta_derivatives(polar_values, cosine_theta_nodes)

Construct `partial_theta Phi` by replacing `Y_lm` with `partial_theta Y_lm`
in the definition of `Phi`. Compute `partial_theta Y_lm` using the recurrence relation
for associated Legendre polynomials. Gauss-Legendre nodes avoid poles.
"""
function build_polar_theta_derivatives(polar_values, cosine_theta_nodes)
    all(x -> -1 < x < 1, cosine_theta_nodes) ||
        throw(ArgumentError("require polar nodes strictly between -1 and 1"))
    sine_theta = sqrt.(1 .- cosine_theta_nodes.^2)

    return Dict((l, m) => begin
        lower_term = l == m ? zero.(values) :
            sqrt((2l + 1) * (l^2 - m^2) / (2l - 1)) .* polar_values[(l-1, m)]
        (l .* cosine_theta_nodes .* values .- lower_term) ./ sine_theta
    end for ((l, m), values) in polar_values)
end

"""
    reconstruct_two_phase_field_derivatives(radial_values, tortoise_slopes,
        tortoise_radius, spin, cloud_frequencies, cloud_m, polar_values,
        grid, phase1_nodes, phase2_nodes)

Calculate and return the combined field and each derivative term inside the
stress-energy tensor of a two-field configuration. Uses `build_polar_theta_derivatives`
to construct `partial_theta Phi` for brevity.
"""
function reconstruct_two_phase_field_derivatives(radial_values, tortoise_slopes,
        tortoise_radius, spin, cloud_frequencies, cloud_m, polar_values,
        grid, phase1_nodes, phase2_nodes)
    length(cloud_frequencies) == length(cloud_m) == 2 ||
        throw(DimensionMismatch("require two carrier frequencies and cloud m values"))
    all(value -> value isa Real && isfinite(value), cloud_frequencies) ||
        throw(ArgumentError("supply finite real carrier frequencies explicitly"))
    keys(radial_values) == keys(tortoise_slopes) ||
        throw(ArgumentError("profiles and slopes must have the same (l, q) keys"))

    geometry = tortoise_geometry(tortoise_radius, spin)
    time_values = Dict((l, q) => -im * sum(q .* cloud_frequencies) * value
        for ((l, q), value) in radial_values)
    azimuthal_values = Dict((l, q) => im * sum(q .* cloud_m) * value
        for ((l, q), value) in radial_values)
    radial_derivatives = Dict(pair => slope / geometry.dr_drstar
        for (pair, slope) in tortoise_slopes)
    polar_theta_derivatives = build_polar_theta_derivatives(polar_values, grid.cosine_theta_nodes)

    reconstruct = (values, polar_table) -> reconstruct_two_phase_field_at_radii(
        values, cloud_m, polar_table, phase1_nodes, phase2_nodes)

    return (; field=reconstruct(radial_values, polar_values),
        derivative_time=reconstruct(time_values, polar_values),
        derivative_radius=reconstruct(radial_derivatives, polar_values),
        derivative_theta=reconstruct(radial_values, polar_theta_derivatives),
        derivative_azimuth=reconstruct(azimuthal_values, polar_values))
end