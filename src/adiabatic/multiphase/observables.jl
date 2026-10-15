# Units: M=1, with the overall Fa^2 divided out of the stress tensor and charges.
# Tensor indices follow (t, r, theta, varphi), even though radial samples use r_*.
# Local radial fluxes are signed currents: positive means toward increasing r.

# Potential-only contribution, retained as a separate check of the full energy.
function two_phase_potential_energy_density(radial_values, tortoise_radius,
        spin, alpha, cloud_m, polar_values, grid, phase1_nodes, phase2_nodes;
        potential=:cosine)
    isempty(phase1_nodes) || isempty(phase2_nodes) ?
        throw(ArgumentError("phase grids must be nonempty")) : nothing

    geometry = tortoise_geometry(tortoise_radius, spin)
    field = reconstruct_two_phase_field_at_radii(radial_values, cloud_m,
        polar_values, phase1_nodes, phase2_nodes)

    # Average U of the joint field, so interactions between components are included.
    # Array axes are (polar node, first phase, second phase).
    phase_count = length(phase1_nodes) * length(phase2_nodes)
    averaged_potential = vec(sum(potential_value.(field, Ref(potential));
        dims=(2, 3))) ./ phase_count

    rho_squared = geometry.radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2
    # 2pi integrates the spatial azimuth; dr/dr_* converts the radial density.
    return 2pi * geometry.dr_drstar * alpha^2 *
        sum(grid.polar_weights .* rho_squared .* averaged_potential)
end

"""
    scalar_stress_tensor(field, derivatives, metric, inverse_metric, alpha;
                         potential=:cosine)

Return T_mu_nu/Fa^2 at one sample point, using signature (-+++).
Supply derivatives in the order [partial_t Phi, partial_r Phi,
partial_theta Phi, partial_varphi Phi].
"""
function scalar_stress_tensor(field::Real, derivatives::AbstractVector{<:Real},
        metric, inverse_metric, alpha; potential=:cosine)
    length(derivatives) == 4 && size(metric) == size(inverse_metric) == (4, 4) ||
        throw(DimensionMismatch("require four derivatives and two 4x4 metrics"))

    # Contract the gradient with the inverse metric, then add V=alpha^2 U(Phi).
    scalar_term = sum(derivatives .* (inverse_metric * derivatives)) / 2 +
        alpha^2 * potential_value(field, potential)

    # T_mu_nu = (partial_mu Phi)(partial_nu Phi) - g_mu_nu * scalar_term.
    return derivatives * transpose(derivatives) - metric * scalar_term
end

"""
    two_phase_observables_at_radius(radial_values, tortoise_slopes,
        tortoise_radius, spin, alpha, cloud_frequencies, cloud_m,
        polar_values, grid, phase1_nodes, phase2_nodes; potential=:cosine)

Return energy and angular momentum per unit tortoise radius, plus the two
signed radial fluxes through that sphere. Integrate over the polar coordinate
and spatial azimuth, and average over both independent phases.

Supply uniform phase grids and real carrier frequencies. The reconstructed
time derivative holds the profiles fixed, omitting slow amplitude changes.
"""
function two_phase_observables_at_radius(radial_values, tortoise_slopes,
        tortoise_radius, spin, alpha, cloud_frequencies, cloud_m,
        polar_values, grid, phase1_nodes, phase2_nodes; potential=:cosine)
    !isempty(grid.polar_weights) && !isempty(phase1_nodes) && !isempty(phase2_nodes) ||
        throw(ArgumentError("quadrature grids must be nonempty"))

    # Input dictionaries are keyed by (l, q). Q=dR/dr_* is converted to dR/dr
    # inside reconstruction, so all derivatives match the Boyer-Lindquist metric.
    samples = reconstruct_two_phase_field_derivatives(radial_values, tortoise_slopes,
        tortoise_radius, spin, cloud_frequencies, cloud_m, polar_values,
        grid, phase1_nodes, phase2_nodes)
    size(samples.field, 1) == length(grid.cosine_theta_nodes) == length(grid.polar_weights) ||
        throw(DimensionMismatch("polar samples, nodes, and weights must agree"))

    geometry = tortoise_geometry(tortoise_radius, spin)
    phase_count = length(phase1_nodes) * length(phase2_nodes)
    energy = angular_momentum = energy_flux = angular_momentum_flux = zero(eltype(samples.field))

    for i in eachindex(grid.cosine_theta_nodes)
        cosine_theta = grid.cosine_theta_nodes[i]
        metric = kerr_metric(tortoise_radius, spin, cosine_theta)
        # The measure is rho^2 d(cos theta) dvarphi. The two phase sums are
        # averages, while the spatial azimuthal integral supplies 2pi.
        weight = 2pi * grid.polar_weights[i] *
            (geometry.radius^2 + spin^2*cosine_theta^2) / phase_count

        for j in eachindex(phase1_nodes), k in eachindex(phase2_nodes)
            derivatives = [samples.derivative_time[i,j,k], samples.derivative_radius[i,j,k],
                samples.derivative_theta[i,j,k], samples.derivative_azimuth[i,j,k]]
            tensor = scalar_stress_tensor(samples.field[i,j,k], derivatives,
                metric.metric, metric.inverse_metric, alpha; potential)
            mixed = metric.inverse_metric * tensor  # Raise the first index: T^mu_nu.

            # Charges use -T^t_t and T^t_varphi; radial currents use
            # -T^r_t and T^r_varphi, respectively.
            energy -= weight * mixed[1, 1]
            angular_momentum += weight * mixed[1, 4]
            energy_flux -= weight * mixed[2, 1]
            angular_momentum_flux += weight * mixed[2, 4]
        end
    end

    # Only charge densities acquire dr/dr_*. Fluxes already refer to a sphere.
    return (; energy_per_tortoise_radius=geometry.dr_drstar * energy,
        angular_momentum_per_tortoise_radius=geometry.dr_drstar * angular_momentum,
        radial_energy_flux=energy_flux,
        radial_angular_momentum_flux=angular_momentum_flux)
end

"""
    ordered_radial_branch(branch, pair_count)

Validate a sequence of profile samples (r*, [R, Q=dR/dr*]) from one radial
integration, and orient its samples toward increasing r_*.

Return `(positions, states, spacing)` with positive spacing, leaving the
input branch unchanged.
"""
function ordered_radial_branch(branch, pair_count)
    # A branch is one radial integration profile series, not an amplitude family.
    # Each state column stores all R coefficients, then all Q=dR/dr_* slopes.
    positions = branch.positions
    length(positions) >= 3 && all(isfinite, positions) ||
        throw(ArgumentError("require at least three finite radial positions"))
    size(branch.states) == (2pair_count, length(positions)) ||
        throw(DimensionMismatch("require one [R; Q] column per radial position"))
    all(isfinite, branch.states) || throw(ArgumentError("require finite profile samples"))

    spacings = diff(positions)
    all(spacings .> 0) || all(spacings .< 0) ||
        throw(ArgumentError("radial positions must be strictly monotonic"))
    # The outer integration normally runs inward. Reverse columns together with
    # their positions so either profile series can be integrated with positive spacing.
    order = first(spacings) > 0 ? collect(eachindex(positions)) :
        reverse(collect(eachindex(positions)))
    ordered_positions = positions[order]
    spacing = (last(ordered_positions) - first(ordered_positions)) / (length(positions)-1)
    all(step -> isapprox(step, spacing; rtol=1e-10, atol=0), diff(ordered_positions)) ||
        throw(ArgumentError("radial quadrature requires uniform branch spacing"))

    return (; positions=ordered_positions, states=branch.states[:, order], spacing)
end

"""
    integrate_two_phase_branch_observables(branch, pairs, spin, alpha,
        cloud_frequencies, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential=:cosine)

Integrate E and J over one profile series prepared by ordered_radial_branch.
Return those charges and the local observables at its lower/upper endpoints.
The (l, q) pairs must follow the ordering of the profile rows.
"""
function integrate_two_phase_branch_observables(branch, pairs, spin, alpha,
        cloud_frequencies, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential=:cosine)
    count = length(pairs)
    # Translate each [R; Q] column into the dictionaries used by reconstruction.
    samples = [two_phase_observables_at_radius(
        Dict(pair => branch.states[index, sample] for (index, pair) in enumerate(pairs)),
        Dict(pair => branch.states[count+index, sample] for (index, pair) in enumerate(pairs)),
        position, spin, alpha, cloud_frequencies, cloud_m, polar_values,
        grid, phase1_nodes, phase2_nodes; potential)
        for (sample, position) in enumerate(branch.positions)]
    all(value -> all(isfinite, values(value)), samples) ||
        throw(ArgumentError("nonfinite observables; check profile magnitudes and quadrature"))

    # Reuse the shared Simpson quadrature. Keep currents at the endpoints;
    # currents are not integrated over radius when computing boundary losses.
    energy = integrate_uniform(
        [value.energy_per_tortoise_radius for value in samples], branch.spacing)
    angular_momentum = integrate_uniform(
        [value.angular_momentum_per_tortoise_radius for value in samples], branch.spacing)
    return (; energy, angular_momentum, lower=first(samples), upper=last(samples))
end

"""
    two_phase_integrated_observables(profiles, spin, alpha,
        cloud_frequencies, components, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential=:cosine)

Combine the two profile series from two_phase_matching_residual(save_profile=true).
Return total E, J over their finite radial domain and signed losses at both
boundaries. The caller is responsible for profile convergence and any tails
outside the retained domain; this function only checks the integration layout.
"""
function two_phase_integrated_observables(profiles, spin, alpha,
        cloud_frequencies, components, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential=:cosine)
    # Match the solver's row order: component first, then retained l within it.
    pairs = [(l, component.q) for component in components for l in component.l_values]
    count = length(pairs)
    count > 0 || throw(ArgumentError("require at least one retained pair"))

    inner = ordered_radial_branch(profiles.from_horizon, count)
    outer = ordered_radial_branch(profiles.from_outer, count)
    # The intervals touch at the match point. Integrating each interval once
    # gives the full domain even though both profile series store that endpoint.
    isapprox(last(inner.positions), first(outer.positions); rtol=0, atol=1e-12) ||
        throw(ArgumentError("radial branches must meet at the same matching position"))

    from_horizon = integrate_two_phase_branch_observables(inner, pairs, spin, alpha,
        cloud_frequencies, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential)
    from_outer = integrate_two_phase_branch_observables(outer, pairs, spin, alpha,
        cloud_frequencies, cloud_m, polar_values, grid,
        phase1_nodes, phase2_nodes; potential)
    # Loss through the inner boundary is toward decreasing r, hence the minus
    # sign. A negative horizon energy loss represents energy gained by the cloud.
    return (; energy=from_horizon.energy + from_outer.energy,
        angular_momentum=from_horizon.angular_momentum + from_outer.angular_momentum,
        horizon_energy_flux=-from_horizon.lower.radial_energy_flux,
        outer_energy_flux=from_outer.upper.radial_energy_flux,
        horizon_angular_momentum_flux=-from_horizon.lower.radial_angular_momentum_flux,
        outer_angular_momentum_flux=from_outer.upper.radial_angular_momentum_flux)
end
