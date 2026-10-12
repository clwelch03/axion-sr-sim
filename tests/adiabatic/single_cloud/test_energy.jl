@testset "Energy density: Schwarzschild without potential" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(1, grid.azimuth_nodes)

    check_position = 10.0
    check_frequency = 0.4
    check_modes = [(1, 1)]
    check_state = ComplexF64[0.7 + 0.2im, -0.03 + 0.01im]

    geometry = tortoise_geometry(check_position, 0.0)
    radial_value, radial_slope = check_state

    expected = 2π * (
        geometry.radius^2 * abs2(check_frequency * radial_value) +
        geometry.radius^2 * abs2(radial_slope) +
        2 * geometry.dr_drstar * abs2(radial_value)
    )

    computed = energy_per_tortoise_radius(
        check_state, check_position, check_frequency,
        0.0, 0.0, check_modes,
        polar_values, azimuthal_values, grid,
    )

    @test isapprox(computed, expected; rtol=1e-12)
end

@testset "Energy density: Kerr without potential" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(1, grid.azimuth_nodes)

    spin = 0.99
    position = 5.0
    frequency = 0.4 + 0.02im
    modes = [(1, 1)]
    state = ComplexF64[0.7 + 0.2im, -0.03 + 0.01im]

    geometry = tortoise_geometry(position, spin)
    radial_factor = geometry.radial_factor
    jacobian = geometry.dr_drstar
    radial_value, radial_slope = state

    expected = 2π * (
        (radial_factor - (4 / 5) * jacobian * spin^2) *
            abs2(frequency * radial_value) +
        radial_factor * abs2(radial_slope) +
        (2jacobian - spin^2 / radial_factor) * abs2(radial_value)
    )

    computed = energy_per_tortoise_radius(
        state, position, frequency, spin, 0.0, modes,
        polar_values, azimuthal_values, grid,
    )

    @test isapprox(computed, expected; rtol=1e-12)
end

@testset "Energy density: finite-amplitude cosine potential" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(1, grid.azimuth_nodes)

    position = 10.0
    alpha = 0.42
    modes = [(1, 1)]
    geometry = tortoise_geometry(position, 0.0)

    for amplitude in (0.2, 0.7, 1.5)
        state = ComplexF64[amplitude, 0.0]

        with_potential = energy_per_tortoise_radius(
            state, position, 0.0, 0.0, alpha, modes,
            polar_values, azimuthal_values, grid,
        )
        without_potential = energy_per_tortoise_radius(
            state, position, 0.0, 0.0, 0.0, modes,
            polar_values, azimuthal_values, grid,
        )

        field_amplitude = sqrt(3) * amplitude
        expected = 4π * geometry.dr_drstar * geometry.radius^2 *
            alpha^2 * (1 - sin(field_amplitude) / field_amplitude)

        @test isapprox(
            with_potential - without_potential, expected; rtol=1e-12,
        )
    end
end

@testset "Energy quadrature: cubic density and branch direction" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(1, grid.azimuth_nodes)
    modes = [(1, 1)]

    start, stop = 2.0, 4.0
    expected = (stop - start) + (stop^4 - start^4) / 4

    for interval_count in (2, 8)
        positions = collect(range(start, stop; length=interval_count + 1))
        states = zeros(ComplexF64, 2, length(positions))

        # Set R=0 and prescribe Q so that the energy density is 1+r_*³.
        # Only the radial term, 2pi*r²*|Q|², contributes.
        for (index, position) in enumerate(positions)
            geometry = tortoise_geometry(position, 0.0)
            states[2, index] = sqrt(
                (1 + position^3) / (2π * geometry.radius^2)
            )
        end

        forward = (; positions, states)
        backward = (
            positions=reverse(positions),
            states=reverse(states; dims=2),
        )

        for branch in (forward, backward)
            computed = integrate_energy_branch(
                branch, 0.0, 0.0, 0.0, modes,
                polar_values, azimuthal_values, grid,
            )

            @test isapprox(computed, expected; rtol=1e-12)
        end
    end
end
