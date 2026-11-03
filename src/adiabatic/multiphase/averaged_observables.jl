# Quadratic stress-tensor terms are averaged analytically. The nonlinear
# potential retains the original finite phase/polar quadrature.

struct PhaseAveragedRadiusWeights
    radial_factor::Float64
    dr_drstar::Float64
    temporal::Vector{Float64}
    azimuthal::Vector{Float64}
    dragging::Float64
    potential_weights::Vector{Float64}
end

"""
Fixed-background data owned by one observable/response calculation. No global
cache or checkpoint fields are introduced. Treat its basis and grids as
read-only; frequencies and radial coefficients remain differentiable inputs.
Phase labels have arbitrary length; the profile solver still has two clouds.
"""
struct PhaseAveragedObservableCache{Q,Basis}
    spin::Float64
    alpha::Float64
    polar_blocks::Vector{Matrix{Float64}}
    theta_blocks::Vector{Matrix{Float64}}
    index_blocks::Vector{UnitRange{Int}}
    phase_vectors::Vector{Q}
    azimuthal_numbers::Vector{Int}
    polar_weights::Vector{Float64}
    cosine_theta::Vector{Float64}
    basis::Basis
    radial_weights::Dict{Float64,PhaseAveragedRadiusWeights}
end

"""
    phase_averaged_observable_cache(spin, alpha, components, cloud_m,
        polar_values, grid, phase_nodes, source_basis)
    phase_averaged_observable_cache(problem::TwoPhaseProblem)

Prepare angular and geometric data for analytic quadratic phase averages.
Supply one node array per independent phase, distinct nonzero/nonconjugate
component labels, and the matching source basis. All angular harmonics sharing
q belong in one component so their interference is retained. Reject phase
aliasing rather than silently changing the original finite-grid model.
"""
function phase_averaged_observable_cache(spin, alpha, components, cloud_m,
        polar_values, grid, phase_nodes, source_basis)
    phase_count = length(cloud_m)
    phase_count > 0 && length(phase_nodes) == phase_count ||
        throw(DimensionMismatch("require one grid per independent phase"))
    !isempty(components) && all(nodes -> !isempty(nodes) && all(isfinite, nodes), phase_nodes) ||
        throw(ArgumentError("require components and nonempty finite phase grids"))
    phase_vectors = [Tuple(component.q) for component in components]
    all(q -> length(q) == phase_count && all(value -> value isa Integer, q), phase_vectors) ||
        throw(DimensionMismatch("require one integer q entry per independent phase"))
    phase_average(q) = prod(sum(cis(-q[i]*node) for node in phase_nodes[i]) /
        length(phase_nodes[i]) for i in 1:phase_count)
    for i in eachindex(phase_vectors), j in eachindex(phase_vectors)
        difference = phase_vectors[i] .- phase_vectors[j]
        combined = phase_vectors[i] .+ phase_vectors[j]
        abs(phase_average(difference) - (i == j ? 1.0 : 0.0)) < 5e-13 &&
            abs(phase_average(combined)) < 5e-13 ||
            throw(ArgumentError("phase grid must resolve quadratic products of distinct, nonzero, nonconjugate components"))
    end
    labels = [(l, component.q) for component in components for l in component.l_values]
    labels == source_basis.mode_labels || throw(ArgumentError("observable/source mode order differs"))
    polar_count = length(grid.polar_weights)
    polar_count > 0 && length(grid.cosine_theta_nodes) == polar_count ||
        throw(DimensionMismatch("polar nodes and weights must match"))
    sample_count = polar_count * prod(length, phase_nodes)
    size(source_basis.basis) == (sample_count, 2length(labels)) &&
        length(source_basis.weights) == length(source_basis.cosine_theta_squared) == sample_count ||
        throw(DimensionMismatch("source basis must match angular/phase grids"))
    theta_table = build_polar_theta_derivatives(polar_values, grid.cosine_theta_nodes)
    polar_blocks, theta_blocks = Matrix{Float64}[], Matrix{Float64}[]
    index_blocks, azimuthal_numbers = UnitRange{Int}[], Int[]
    offset = 0
    for component in components
        m = sum(component.q .* cloud_m)
        !isempty(component.l_values) || throw(ArgumentError("require angular modes in each component"))
        polar = hcat([polar_values[(l, abs(m))] for l in component.l_values]...)
        size(polar, 1) == polar_count || throw(DimensionMismatch("polar basis must match quadrature"))
        push!(polar_blocks, polar)
        push!(theta_blocks, hcat([theta_table[(l, abs(m))] for l in component.l_values]...))
        push!(index_blocks, offset+1:offset+length(component.l_values))
        push!(azimuthal_numbers, m)
        offset += length(component.l_values)
    end
    return PhaseAveragedObservableCache(spin, alpha, polar_blocks, theta_blocks,
        index_blocks, phase_vectors, azimuthal_numbers, collect(grid.polar_weights),
        collect(grid.cosine_theta_nodes), source_basis, Dict{Float64,PhaseAveragedRadiusWeights}())
end

phase_averaged_observable_cache(problem::TwoPhaseProblem) = phase_averaged_observable_cache(
    problem.spin, problem.alpha, problem.components, problem.cloud_m,
    problem.polar_values, problem.grid, (problem.phase1_nodes, problem.phase2_nodes), problem.source_basis)

function phase_averaged_radius_weights(position, cache::PhaseAveragedObservableCache)
    return get!(cache.radial_weights, Float64(position)) do
        geometry = tortoise_geometry(position, cache.spin)
        sine_squared = 1 .- cache.cosine_theta.^2
        azimuth_factor = geometry.radial_factor^2 .- cache.spin^2 * geometry.delta .* sine_squared
        temporal = azimuth_factor ./ geometry.radial_factor
        azimuthal = (geometry.delta .- cache.spin^2 .* sine_squared) ./
            (geometry.radial_factor .* sine_squared)
        dragging = 2cache.spin * geometry.radius / geometry.radial_factor
        # basis.weights already includes the polar weight and 1/(n1*n2).
        potential_weights = (2pi * cache.alpha^2 * geometry.dr_drstar) .*
            cache.basis.weights .* (geometry.radius^2 .+
                cache.spin^2 .* cache.basis.cosine_theta_squared)
        PhaseAveragedRadiusWeights(geometry.radial_factor, geometry.dr_drstar,
            temporal, azimuthal, dragging, potential_weights)
    end
end

function phase_averaged_potential_energy(coefficients, basis, weights, potential)
    field = 2 .* (basis.basis * coefficients)
    return sum(weights .* potential_value.(field, Ref(potential)))
end

# Same discrete potential quadrature, with an explicit local chain rule to
# avoid dense Float64-by-Dual multiplication at every radius. U'=Phi-h(Phi).
# Nested Duals and other precisions use the generic method above.
function phase_averaged_potential_energy(
        coefficients::AbstractVector{ForwardDiff.Dual{Tag,Float64,N}},
        basis, weights, potential) where {Tag,N}
    values = ForwardDiff.value.(coefficients)
    field = 2 .* (basis.basis * values)
    energy = sum(weights .* potential_value.(field, Ref(potential)))
    potential_derivative = field .- potential_remainder.(field, Ref(potential))
    gradient = 2 .* (transpose(basis.basis) * (weights .* potential_derivative))
    directions = ntuple(direction -> sum(gradient[index] *
        ForwardDiff.partials(coefficients[index])[direction] for index in eachindex(coefficients)), Val(N))
    return ForwardDiff.Dual{Tag}(energy, directions)
end

"""
Evaluate phase-averaged densities and currents on the original polar/phase grid.

For each q, form the full polar profiles F_q, Q_q=dF_q/dr*, and dF_q/dtheta
before taking products; this preserves all same-q cross terms between l values.
Resolved phase orthogonality gives <Da Phi Db Phi>=2Re sum_q Da,q conj(Db,q).

With S=r^2+a^2 and B=S^2-a^2 Delta sin^2(theta), the quadratic energy density
per dr* is 2pi sum_i w_i[(B/S omega_q^2 + (Delta-a^2 sin^2)/(S sin^2) m_q^2)
|F_q|^2 + S|Q_q|^2 + (Delta/S)|dtheta F_q|^2]. The nonlinear potential is
evaluated on exactly the original finite quadrature. This follows by canceling
the time-azimuth cross term in -T^t_t before evaluation, not by dropping it.
"""
function phase_averaged_observables_at_radius(radial_values, tortoise_slopes,
        position, carrier_frequencies, cache::PhaseAveragedObservableCache; potential=:cosine)
    length(carrier_frequencies) == length(first(cache.phase_vectors)) ||
        throw(DimensionMismatch("require one carrier frequency per independent phase"))
    all(value -> value isa Real && isfinite(value), carrier_frequencies) ||
        throw(ArgumentError("supply finite real carrier frequencies"))
    length(radial_values) == length(tortoise_slopes) == length(cache.basis.mode_labels) ||
        throw(DimensionMismatch("profile and slope ordering must match the source basis"))
    weights = phase_averaged_radius_weights(position, cache)
    scalar_zero = zero(real(first(radial_values)) + first(carrier_frequencies))
    energy = angular_momentum = energy_flux = angular_momentum_flux = scalar_zero
    for component_index in eachindex(cache.index_blocks)
        indices = cache.index_blocks[component_index]
        polar = cache.polar_blocks[component_index]
        theta = cache.theta_blocks[component_index]
        frequency = sum(cache.phase_vectors[component_index] .* carrier_frequencies)
        azimuthal_number = cache.azimuthal_numbers[component_index]
        for polar_index in eachindex(cache.polar_weights)
            field = slope = theta_derivative = zero(first(radial_values))
            for (column, radial_index) in enumerate(indices)
                field += polar[polar_index, column] * radial_values[radial_index]
                slope += polar[polar_index, column] * tortoise_slopes[radial_index]
                theta_derivative += theta[polar_index, column] * radial_values[radial_index]
            end
            field_squared = abs2(field)
            polar_weight = cache.polar_weights[polar_index]
            energy += (2pi * polar_weight) *
                ((weights.temporal[polar_index] * frequency^2 +
                  weights.azimuthal[polar_index] * azimuthal_number^2) * field_squared +
                 weights.radial_factor * abs2(slope) + weights.dr_drstar * abs2(theta_derivative))
            angular_momentum += (4pi * polar_weight) *
                (weights.temporal[polar_index] * frequency * azimuthal_number -
                 weights.dragging * azimuthal_number^2) * field_squared
            current = (4pi * weights.radial_factor * polar_weight) * imag(conj(field) * slope)
            energy_flux += frequency * current
            angular_momentum_flux += azimuthal_number * current
        end
    end
    coefficients = vcat(real.(radial_values), imag.(radial_values))
    energy += phase_averaged_potential_energy(coefficients, cache.basis,
        weights.potential_weights, potential)
    return (; energy_per_tortoise_radius=energy,
        angular_momentum_per_tortoise_radius=angular_momentum,
        radial_energy_flux=energy_flux, radial_angular_momentum_flux=angular_momentum_flux)
end

"""
    phase_averaged_integrated_observables(radial, carrier_frequencies, cache;
        potential=:cosine)

Integrate the averaged densities using the existing radial quadrature and
boundary-loss signs. The caller remains responsible for profile convergence.
"""
function phase_averaged_integrated_observables(radial, carrier_frequencies,
        cache::PhaseAveragedObservableCache; potential=:cosine)
    count = length(cache.basis.mode_labels)
    inner = ordered_radial_branch(radial.from_horizon, count)
    outer = ordered_radial_branch(radial.from_outer, count)
    isapprox(last(inner.positions), first(outer.positions); rtol=0, atol=1e-12) ||
        throw(ArgumentError("radial branches must meet at the same matching position"))
    function integrate_branch(branch)
        samples = [phase_averaged_observables_at_radius(
            view(branch.states, 1:count, i), view(branch.states, count+1:2count, i),
            position, carrier_frequencies, cache; potential)
            for (i, position) in enumerate(branch.positions)]
        all(value -> all(isfinite, values(value)), samples) ||
            throw(ArgumentError("nonfinite observables; check profile magnitudes and quadrature"))
        energy = integrate_uniform([value.energy_per_tortoise_radius for value in samples], branch.spacing)
        angular_momentum = integrate_uniform([value.angular_momentum_per_tortoise_radius for value in samples], branch.spacing)
        return (; energy, angular_momentum, lower=first(samples), upper=last(samples))
    end
    horizon, infinity = integrate_branch(inner), integrate_branch(outer)
    return (; energy=horizon.energy+infinity.energy,
        angular_momentum=horizon.angular_momentum+infinity.angular_momentum,
        horizon_energy_flux=-horizon.lower.radial_energy_flux,
        outer_energy_flux=infinity.upper.radial_energy_flux,
        horizon_angular_momentum_flux=-horizon.lower.radial_angular_momentum_flux,
        outer_angular_momentum_flux=infinity.upper.radial_angular_momentum_flux)
end
