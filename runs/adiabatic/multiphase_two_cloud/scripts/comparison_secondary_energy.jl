# Comparison coordinate: quadratic secondary energy on the averaged primary.
# This extends the cosine diagnostic to other potentials without changing the
# faithful Omiya solver or assigning mixed-mode energy to either bound cloud.
function comparison_secondary_energy_density(state, position, primary_values,
        frequency, spin, alpha, modes, polar, azimuthal, grid; potential=:cosine)
    if potential === :cosine
        return secondary_energy_per_tortoise_radius(state, position, primary_values,
            frequency, spin, alpha, modes, polar, azimuthal, grid)
    end
    kinetic = energy_per_tortoise_radius(state, position, frequency, spin, 0.0,
        modes, polar, azimuthal, grid; fundamental_m=last(first(modes)))
    field = reconstruct_field(primary_values, polar, azimuthal, grid)
    # h(phi)=phi-U'(phi), hence U''=1-h'. This differentiates the same stable
    # potential remainder used by the nonlinear radial source.
    curvature = map(value -> 1 - ForwardDiff.derivative(
        x -> potential_remainder(x, potential), value), field)
    averaged_curvature = azimuthal_average(curvature, grid)
    geometry = tortoise_geometry(position, spin)
    rho_squared = geometry.radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2
    profile = sum(state[i] .* polar[mode] for (i, mode) in enumerate(modes))
    mass_energy = 2pi * geometry.dr_drstar * alpha^2 * sum(
        grid.polar_weights .* rho_squared .* averaged_curvature .* abs2.(profile))
    return kinetic + mass_energy
end

comparison_energy_convention(potential) =
    "isolated $(potential) primary; quadratic secondary using averaged U'' of $(potential); mixed modes excluded from labels"
