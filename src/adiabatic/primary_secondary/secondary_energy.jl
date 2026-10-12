"""
    secondary_energy_per_tortoise_radius(state, tortoise_radius, primary_values,
                                         frequency, spin, alpha, secondary_modes,
                                         polar_values, azimuthal_values, grid)

Return the averaged secondary's energy per unit tortoise radius, `d(E2/Fa²)/dr_*`.

The secondary's potential energy is the second-order term in ϕ₂ of
`mu^2 (1 - cos(ϕ₁ + ϕ₂))`, period-averaged: `(mu^2/2) <cos ϕ₁> ϕ₂²`, a mass
term with effective mass squared `mu^2 <cos ϕ₁>`, consistent with the
averaged Eq. (61). Kinetic and gradient terms come from
[`energy_per_tortoise_radius`](@ref).
"""
function secondary_energy_per_tortoise_radius(state, tortoise_radius, primary_values,
                                              frequency, spin, alpha, secondary_modes,
                                              polar_values, azimuthal_values, grid)
    secondary_m = last(first(secondary_modes))

    # alpha enters energy_per_tortoise_radius only through its potential term
    # alpha^2 U(ϕ), so alpha = 0 leaves the kinetic and gradient energy.
    kinetic_and_gradient = energy_per_tortoise_radius(
        state, tortoise_radius, frequency, spin, 0.0, secondary_modes,
        polar_values, azimuthal_values, grid; fundamental_m=secondary_m,
    )

    geometry = tortoise_geometry(tortoise_radius, spin)
    rho_squared = geometry.radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2
    mass_deficit_average = averaged_mass_deficit(
        primary_values, polar_values, azimuthal_values, grid,
    )
    # The secondary's polar profile sum_l R_l Y_l(x); its azimuthal average of
    # ϕ₂² is 2|profile|², and the azimuthal integral contributes 2π.
    secondary_profile = sum(
        state[index] .* polar_values[mode] for (index, mode) in enumerate(secondary_modes)
    )
    # (alpha^2/2) <cos ϕ₁> <ϕ₂²>, with <cos ϕ₁> = 1 - <1 - cos ϕ₁>.
    mass_energy = 2π * geometry.dr_drstar * alpha^2 * sum(
        grid.polar_weights .* rho_squared .* (1 .- mass_deficit_average) .*
        abs2.(secondary_profile)
    )

    return kinetic_and_gradient + mass_energy
end

# Outward energy current through a sphere at one tortoise radius, for modes
# sharing the single frequency omega2 (the conjugate partners included).
function outward_secondary_flux(state, tortoise_radius, frequency, spin, mode_count)
    geometry = tortoise_geometry(tortoise_radius, spin)
    return 4π * geometry.radial_factor * sum(
        real(im * frequency * state[index] * conj(state[mode_count + index]))
        for index in 1:mode_count
    )
end

"""
    secondary_energy_and_fluxes(solution, primary_profiles, primary_modes, spin, alpha,
                                secondary_modes, couplings, polar_values,
                                azimuthal_values, grid)

Return `(; energy, horizon_flux, outer_flux)` for an averaged secondary
`solution` from [`solve_secondary_mode`](@ref), in units `Fa² M`.
Positive flux means energy lost by the secondary. For a converged solution,
`2 imag(omega2) * energy ≈ -(horizon_flux + outer_flux)`.
"""
function secondary_energy_and_fluxes(solution, primary_profiles, primary_modes, spin, alpha,
                                     secondary_modes, couplings, polar_values,
                                     azimuthal_values, grid)
    secondary = secondary_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        primary_profiles, primary_modes, spin, alpha, secondary_modes, couplings,
        polar_values, azimuthal_values, grid; save_profile=true, reading=:averaged,
    )

    energy = 0.0
    branch_pairs = ((secondary.from_horizon, primary_profiles.from_horizon),
                    (secondary.from_outer, primary_profiles.from_outer))
    for (secondary_part, primary_part) in branch_pairs
        positions = secondary_part.positions
        # Integrate toward increasing tortoise radius.
        order = issorted(positions) ? eachindex(positions) : reverse(eachindex(positions))
        densities = [
            secondary_energy_per_tortoise_radius(
                secondary_part.states[:, sample], positions[sample],
                primary_values_at(primary_part, positions[sample], primary_modes),
                solution.frequency, spin, alpha, secondary_modes,
                polar_values, azimuthal_values, grid,
            )
            for sample in order
        ]
        energy += integrate_uniform(densities, abs(positions[2] - positions[1]))
    end

    mode_count = length(secondary_modes)
    # Loss into the horizon is minus the outward current at the inner endpoint.
    horizon_flux = -outward_secondary_flux(
        secondary.from_horizon.states[:, 1], first(secondary.from_horizon.positions),
        solution.frequency, spin, mode_count,
    )
    outer_flux = outward_secondary_flux(
        secondary.from_outer.states[:, 1], first(secondary.from_outer.positions),
        solution.frequency, spin, mode_count,
    )
    return (; energy, horizon_flux, outer_flux)
end
