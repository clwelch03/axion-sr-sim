@testset "Polynomial potentials: analytic source and energy" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(1, grid.azimuth_nodes)

    modes = [(1, 1)]
    amplitude = 1.3 + 0.4im
    radial_values = Dict((1, 1) => amplitude)
    state = ComplexF64[amplitude, 0.0]
    position, alpha = 10.0, 0.42

    for spin in (0.0, 0.99)
        geometry = tortoise_geometry(position, spin)
        radius_squared = geometry.radius^2

        # Exact polar moments: ∫Y11⁴ dx = 3/5,
        # ∫Y11⁶ dx = 27/70.
        cubic_source = alpha^2 * abs2(amplitude) * amplitude *
            (3radius_squared / 10 + 3spin^2 / 70)

        quintic_source = alpha^2 * abs2(amplitude)^2 * amplitude *
            (9 / 280) * (radius_squared + spin^2 / 9)

        # For a pure Y11 field, b = sqrt(3)*|A| and
        # ∫phi^(2n) dΩ = 4π*b^(2n)/(2n+1).
        # An extra cos²θ multiplies this moment by 1/(2n+3).
        field_amplitude = sqrt(3) * abs(amplitude)

        quartic_integral =
            field_amplitude^2 / 6 *
                (radius_squared + spin^2 / 5) -
            field_amplitude^4 / 120 *
                (radius_squared + spin^2 / 7)

        sextic_addition = field_amplitude^6 / 5040 *
            (radius_squared + spin^2 / 9)

        for potential in (:quartic, :sextic)
            source = nonlinear_source(
                radial_values, geometry.radius, spin, alpha,
                polar, azimuthal, grid;
                modes, potential,
            )

            expected_source = potential === :quartic ?
                cubic_source : cubic_source - quintic_source

            @test isapprox(
                source[(1, 1)], expected_source; rtol=1e-12,
            )

            with_potential = energy_per_tortoise_radius(
                state, position, 0.4, spin, alpha, modes,
                polar, azimuthal, grid; potential,
            )
            without_potential = energy_per_tortoise_radius(
                state, position, 0.4, spin, 0.0, modes,
                polar, azimuthal, grid; potential,
            )

            angular_integral = potential === :quartic ?
                quartic_integral : quartic_integral + sextic_addition

            expected_energy = 4π * geometry.dr_drstar *
                alpha^2 * angular_integral

            @test isapprox(
                with_potential - without_potential,
                expected_energy; rtol=1e-12,
            )
        end
    end
end