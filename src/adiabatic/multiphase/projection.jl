"""
    reconstruct_two_phase_field_at_radius(radial_values, cloud_m, polar_values,
                                        phase1_nodes, phase2_nodes)

Reconstruct the real field at one radius for two independent cloud phases.

Supply complex radial coefficients keyed by (l, q), with nonzero integer
2-tuples q and one representative per conjugate pair. `cloud_m` lists the two
cloud azimuthal integers.

Return a real array with axes (polar node, first phase, second phase), using
2*real(sum(R_lq * Y_l,abs(m_q) * exp(-im*q⋅psi))).
"""
function reconstruct_two_phase_field_at_radius(radial_field_values, cloud_m_vec, polar_field_values,
                                     phase1_nodes, phase2_nodes)
    length(cloud_m_vec) == 2 || throw(DimensionMismatch("require two cloud m values"))
    field_type = promote_type(Float64,
        (typeof(real(value)) for value in values(radial_field_values))...)
    field = zeros(field_type, length(first(values(polar_field_values))),
        length(phase1_nodes), length(phase2_nodes))

    for ((l, q), radial_value) in radial_field_values
        length(q) == 2 || throw(DimensionMismatch("require two entries in q"))
        m = sum(q .* cloud_m_vec)
        polar = polar_field_values[(l, abs(m))]
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
    two_phase_field_component_coefficients(field, output_pairs, cloud_m, polar_values,
                            polar_weights, phase1_nodes, phase2_nodes)

Extract coefficients of requested `(l, q)` angular/phase patterns at one radius.

`field` has axes (polar node, first phase, second phase). Supply integer
2-tuples `q` and the same cloud ordering and polar table as for reconstruction.

Project onto each requested angular/phase basis function by summing field * Y_l,
abs(m_q) * exp(+im*q⋅psi) over all grid nodes, weighted by the polar quadrature
weights and divided by the number of phase pairs. Return the resulting coefficients
in a dictionary keyed by the requested `(l, q)` pairs without an extra factor of 2.
"""
function two_phase_field_component_coefficients(field, output_pairs, cloud_m, polar_values,
                                 polar_weights, phase1_nodes, phase2_nodes)
    # Check input validity
    length(cloud_m) == 2 || throw(DimensionMismatch("require two cloud m values"))
    !isempty(polar_weights) && !isempty(phase1_nodes) && !isempty(phase2_nodes) ||
        throw(ArgumentError("quadrature grids must be nonempty"))
    size(field) == (length(polar_weights), length(phase1_nodes), length(phase2_nodes)) ||
        throw(DimensionMismatch("field axes must match the quadrature grids"))

    # Total number of phase pairs we average over
    phase_count = length(phase1_nodes) * length(phase2_nodes)

    # Calculate the coefficient for one (l, q) pair
    function component_coefficient(l, q)
        length(q) == 2 ||
            throw(DimensionMismatch("require two entries in q"))

        m = sum(q .* cloud_m)
        polar = polar_values[(l, abs(m))]
        length(polar) == length(polar_weights) ||
            throw(DimensionMismatch("polar table and weights must match"))

        sum(polar_weights[i] * polar[i] * field[i, j, k] *
            cis(q[1]*phase1_nodes[j] + q[2]*phase2_nodes[k])
            for i in eachindex(polar_weights), j in eachindex(phase1_nodes),
                k in eachindex(phase2_nodes)) / phase_count
    end

    return Dict((l, q) => component_coefficient(l, q) for (l, q) in output_pairs)
end

"""
    TwoPhaseSourceBasis(mode_labels, basis, weights, cosine_theta_squared)

Fixed angular/phase data reused when evaluating the joint nonlinear source.

- `mode_labels`: the ordered `(l, q)` mode labels.
- `basis`: reconstruction/projection matrix, with cosine columns followed by
  sine columns and rows in flattened polar/phase-sample order.
- `weights`: polar quadrature weights including the two-phase averaging factor.
- `cosine_theta_squared`: squared polar nodes in the same sample order.

Construct a consistent basis with [`build_two_phase_source_basis`](@ref).
The field types are inferred from the supplied data. This container is
immutable, but its arrays remain mutable; treat them as read-only after
construction. The constructor stores the supplied arrays without copying or
validating them. Radial coefficients, field samples, and source values are
computed separately and are not stored here.
"""
struct TwoPhaseSourceBasis{ModeLabels, Basis, Weights, CosineThetaSquared}
    mode_labels::ModeLabels
    basis::Basis
    weights::Weights
    cosine_theta_squared::CosineThetaSquared
end

"""
    build_two_phase_source_basis(mode_labels, cloud_m, polar_values, grid,
                                 phase1_nodes, phase2_nodes)

Precompute the angular/phase basis and quadrature weights for repeated source evaluations.

Supply `mode_labels` as an ordered collection of distinct `(l, q)` labels, with
nonzero integer 2-tuples `q` and one representative per conjugate pair.
`cloud_m` contains the two clouds' azimuthal integers. The real table
`polar_values[(l, abs(m))]`,
where `m = sum(q .* cloud_m)`, must be evaluated at `grid.cosine_theta_nodes`.
Both phase-node arrays must sample one full period uniformly; the polar and
phase grids must be nonempty. These input conventions are assumed, not validated.

For `n = length(mode_labels)` and `p` angular/phase samples, return a
[`TwoPhaseSourceBasis`](@ref) with these fields:

- `mode_labels`: the supplied mode order, collected into a vector.
- `basis`: a `p x 2n` Float64 matrix. For each mode, the first half of columns
  contains its polar harmonic times `cos(q[1]*phase1 + q[2]*phase2)`; the
  second half contains the corresponding sine pattern.
- `weights`: polar quadrature weights repeated for every phase pair and
  divided by `length(phase1_nodes) * length(phase2_nodes)` to average over phases.
- `cosine_theta_squared`: the squared polar nodes repeated for every phase pair.

Rows follow `vec(field)` ordering: polar index first, then the first phase,
then the second phase. Thus `2 * basis * vcat(real.(R), imag.(R))` reconstructs
the flattened real field when the complex coefficients `R` follow `mode_labels`.

Reuse the basis while the mode order, cloud m values, polar table, and grids
remain fixed; rebuild it when those inputs change. It stores fixed basis data,
not a field or source value, and contains no radial coefficients or background
parameters. Treat its arrays as read-only during evaluation.

See also [`nonlinear_two_phase_source`](@ref).
"""
function build_two_phase_source_basis(mode_labels, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes)
    mode_labels = collect(mode_labels)
    polar_node_count = length(grid.polar_weights)
    phase_pair_count = length(phase1_nodes) * length(phase2_nodes)
    num_radial_modes = length(mode_labels)
    basis = Matrix{Float64}(undef, polar_node_count * phase_pair_count, 2num_radial_modes)

    # Consider a single mode with phase combination `phase`
    # phase = q[1] * phase1_nodes[j] + q[2] * phase2_nodes[k]
    # That mode contributes to the field as
    # 2 * (polar harmonic) * [Re(radial coeff) * cos(phase) + Im(radial coeff) * sin(phase)]
    # So we need to store (polar harmonic) * sin(phase) and (polar harmonic) * cos(phase)
    # Each row in `basis` corresponds to a single polar/phase sampling point
    # First (second) half of the columns contains the cosine (sine) patterns
    for (column, (l, q)) in enumerate(mode_labels)
        polar_harmonic_values = polar_values[(l, abs(sum(q .* cloud_m)))]
        for k in eachindex(phase2_nodes), j in eachindex(phase1_nodes)
            phase = q[1] * phase1_nodes[j] + q[2] * phase2_nodes[k]
            for i in eachindex(polar_harmonic_values)
                row = i + polar_node_count * ((j - 1) + length(phase1_nodes) * (k - 1))
                basis[row, column] = polar_harmonic_values[i] * cos(phase)
                basis[row, num_radial_modes + column] = polar_harmonic_values[i] * sin(phase)
            end
        end
    end

    # Use the same flattened sample ordering as the basis rows.
    weights = repeat(grid.polar_weights, phase_pair_count) ./ phase_pair_count
    cosine_theta_squared = repeat(grid.cosine_theta_nodes.^2, phase_pair_count)
    return TwoPhaseSourceBasis(mode_labels, basis, weights, cosine_theta_squared)
end

"""
    nonlinear_two_phase_source(parameters, source_basis; potential=:cosine)

Evaluate the nonlinear source at one radius by reconstruction and projection
with a precomputed basis, accepting and returning packed real vectors.
A prepared basis is required for every evaluation.

Let `n = length(source_basis.mode_labels)`. Supply a real parameter vector of length
`2n + 3`, ordered as `[real.(R); imag.(R); radius; spin; alpha]`, where the
complex radial coefficients `R` follow the basis's mode order. `radius` is
Boyer–Lindquist r in M=1 units. Use a basis from
[`build_two_phase_source_basis`](@ref) and potential `:cosine`, `:quartic`,
or `:sextic`.

Reconstruct the flattened joint field as `2 * basis * coefficients`. Apply
`potential_remainder` to that field, multiply by the precomputed quadrature
weights and `radius^2 + spin^2*cos(theta)^2`, then project with `transpose(basis)`.
Reconstruction includes the conjugate modes through its factor of two;
projection introduces no additional factor of two.

Return a real vector of length `2n`: all real source coefficients followed by
all imaginary source coefficients, in basis order. The result includes
`alpha^2`; the radial RHS applies its minus sign and tortoise factor separately.

This generic method propagates ForwardDiff values through the matrix arithmetic
when they are supplied in `parameters`. It does not explicitly construct the
local source Jacobian. The fixed basis is reused without modification.

See also [`build_two_phase_source_basis`](@ref), [`potential_remainder`](@ref).
"""
function nonlinear_two_phase_source(parameters, source_basis::TwoPhaseSourceBasis; potential=:cosine)
    num_radial_modes = length(source_basis.mode_labels)
    length(parameters) == 2num_radial_modes + 3 ||
        throw(DimensionMismatch("expected real and imaginary coefficients, radius, spin, and alpha"))
    coefficients = parameters[1:2num_radial_modes]
    radius, spin, alpha = parameters[2num_radial_modes+1:2num_radial_modes+3]

    # Reconstruct the real field at every angular/phase sample.
    field = 2 .* (source_basis.basis * coefficients)

    # Apply the potential, geometry, and quadrature weights pointwise.
    geometry_factor = radius^2 .+ spin^2 .* source_basis.cosine_theta_squared
    remainder = potential_remainder.(field, Ref(potential))
    weighted_remainder = source_basis.weights .* geometry_factor .* remainder

    # Return all real source coefficients, followed by all imaginary ones.
    return alpha^2 .* (transpose(source_basis.basis) * weighted_remainder)
end

"""
    local_source_jacobian(parameters, source_basis; potential=:cosine)

Return `(source, jacobian)` for the packed source at one radius,
where the jacobian is `J_{ij} = partial source[i]/partial parameters[j]`.
Supply ordinary real parameters in the order used by `nonlinear_two_phase_source`:
`[real.(R); imag.(R); radius; spin; alpha]`, with R in basis order.
For n retained modes, the source has length 2n and its Jacobian is
2n by (2n + 3); columns follow the parameter order. Hold the basis fixed.
Differentiate the implemented potential remainder, including its small-field
branch. The returned derivatives do not include radial integration or matching.
"""
function local_source_jacobian(parameters, source_basis::TwoPhaseSourceBasis; potential=:cosine)
    mode_count = length(source_basis.mode_labels)
    coefficients = parameters[1:2mode_count]
    radius, spin, alpha = parameters[2mode_count+1:2mode_count+3]
    field = 2 .* (source_basis.basis * coefficients)

    # Evaluate h(field) and its derivative with respect to field together.
    remainder_duals = [
        potential_remainder(ForwardDiff.Dual{Nothing}(value, (1.0,)), potential)
        for value in field
    ]
    remainder = ForwardDiff.value.(remainder_duals)
    remainder_derivative = [ForwardDiff.partials(value)[1] for value in remainder_duals]

    # The source value uses the same quadrature as nonlinear_two_phase_source.
    geometry_weights = source_basis.weights .*
        (radius^2 .+ spin^2 .* source_basis.cosine_theta_squared)
    projected = transpose(source_basis.basis) * (geometry_weights .* remainder)

    # Chain rule: coefficient changes -> field changes -> source changes.
    coefficient_jacobian = (2alpha^2) .* (transpose(source_basis.basis) *
        ((geometry_weights .* remainder_derivative) .* source_basis.basis))

    # Final three columns differentiate radius, spin, and alpha.
    background_columns = hcat(
        2radius .* source_basis.weights .* remainder,
        2spin .* source_basis.weights .* source_basis.cosine_theta_squared .* remainder,
    )
    background_jacobian = alpha^2 .* (transpose(source_basis.basis) * background_columns)
    source = alpha^2 .* projected
    jacobian = hcat(coefficient_jacobian, background_jacobian, 2alpha .* projected)
    return source, jacobian
end

"""
"""
function nonlinear_two_phase_source(
        parameters::AbstractVector{ForwardDiff.Dual{Tag,Float64,N}},
        source_basis::TwoPhaseSourceBasis; potential=:cosine) where {Tag,N}
    length(parameters) == 2length(source_basis.mode_labels) + 3 ||
        throw(DimensionMismatch("expected real and imaginary coefficients, radius, spin, and alpha"))

    # Evaluate the source and its local Jacobian using ordinary real values.
    source, source_jacobian = local_source_jacobian(
        ForwardDiff.value.(parameters), source_basis; potential)

    # Carry the incoming parameter derivatives through the local source map.
    parameter_derivatives = [
        ForwardDiff.partials(parameter)[direction]
        for parameter in parameters, direction in 1:N
    ]
    source_derivatives = source_jacobian * parameter_derivatives

    # Return Dual values with the original differentiation tag and directions.
    return [ForwardDiff.Dual{Tag}(source[index],
        ntuple(direction -> source_derivatives[index, direction], Val(N)))
        for index in eachindex(source)]
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

    reconstruct = (values, polar_table) -> reconstruct_two_phase_field_at_radius(
        values, cloud_m, polar_table, phase1_nodes, phase2_nodes)

    return (; field=reconstruct(radial_values, polar_values),
        derivative_time=reconstruct(time_values, polar_values),
        derivative_radius=reconstruct(radial_derivatives, polar_values),
        derivative_theta=reconstruct(radial_values, polar_theta_derivatives),
        derivative_azimuth=reconstruct(azimuthal_values, polar_values))
end
