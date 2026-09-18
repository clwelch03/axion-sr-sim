@testset "Tortoise geometry" begin
    for radius in (2.1, 3.0, 20.0)
        coordinate = radius + 2log(radius - 2)
        @test tortoise_coordinate(radius - 2, 0.0) ≈ coordinate
        @test tortoise_geometry(coordinate, 0.0).radius ≈ radius
    end
    for spin in (0.0, 0.5, 0.99), coordinate in (-100.0, 5.0, 100.0)
        geometry = tortoise_geometry(coordinate, spin)
        @test isapprox(tortoise_coordinate(geometry.horizon_distance, spin), coordinate; atol=1e-11)
        step = 1e-4
        distance_plus = tortoise_geometry(coordinate + step, spin).horizon_distance
        distance_minus = tortoise_geometry(coordinate - step, spin).horizon_distance
        @test isapprox((distance_plus - distance_minus) / (2step), geometry.dr_drstar; rtol=1e-7)
    end
    deep = tortoise_geometry(-100.0, 0.0)
    @test deep.radius == deep.r_plus  # s is below the spacing of Float64 at r_+.
    @test deep.horizon_distance > 0 && deep.delta > 0
    @test isapprox(deep.log_horizon_distance, -51; atol=1e-12)
    @test_throws DomainError kerr_horizons(1.0)
end

@testset "Tortoise equation from the chain rule" begin
    spin, alpha, frequency = 0.5, 0.42, 0.4088 + 1.5e-7im
    value, radial_slope, mixed_value, source = 0.3 + 0.2im, -0.1 + 0.4im, 0.07im, 0.2 - 0.1im
    for coordinate in (-3.0, 5.0, 20.0)
        geometry = tortoise_geometry(coordinate, spin)
        radial = radial_operator_coefficients(geometry.radius, spin, alpha, frequency, 1, 1)
        tortoise = tortoise_operator_coefficients(geometry, spin, alpha, frequency, 1, 1)
        radial_curvature = radial_second_derivative(value, radial_slope, mixed_value, source, radial)
        h = geometry.dr_drstar
        h_prime = (radial.delta_prime * geometry.radial_factor - 2geometry.radius * radial.delta) /
                  geometry.radial_factor^2
        expected = h^2 * radial_curvature + h * h_prime * radial_slope
        actual = -tortoise.slope_coefficient * h * radial_slope - tortoise.local_coefficient * value -
                 tortoise.coupling_scale * mixed_value - tortoise.source_scale * source
        @test isapprox(actual, expected; rtol=1e-12)
    end
    geometry = tortoise_geometry(-400.0, 0.99)
    coefficients = tortoise_operator_coefficients(geometry, 0.99, alpha, frequency, 1, 1)
    horizon_frequency = 0.99 / (2geometry.r_plus)
    @test coefficients.local_coefficient ≈ (frequency - horizon_frequency)^2
    @test abs(coefficients.coupling_scale) < 1e-20
    @test all(isfinite, values(coefficients))
end

@testset "Coupled first-order state and harmonic frequencies" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    # Even m and an intentionally mixed ordering check m=n*m0, not m=n.
    modes = [(4, 4), (2, 2), (4, 2)]
    state = ComplexF64[0.1im, 0.2 + 0.1im, -0.3, 0.04, 0.05im, -0.06]
    saved_state = copy(state)
    derivative = similar(state)
    linear_tortoise_rhs!(derivative, state, 5.0, 0.5, 0.42, 0.3 + 0.01im,
                         modes, couplings; fundamental_m=2)
    @test state == saved_state
    @test derivative[1:3] == state[4:6]
    geometry = tortoise_geometry(5.0, 0.5)
    for (index, (l, m)) in enumerate(modes)
        radial = radial_operator_coefficients(geometry.radius, 0.5, 0.42, (m / 2) * (0.3 + 0.01im), l, m)
        mixed = sum(couplings[(l, other_l, m)] * state[j]
                    for (j, (other_l, other_m)) in enumerate(modes) if m == other_m)
        h = geometry.dr_drstar
        h_prime = (radial.delta_prime - 2geometry.radius * h) / geometry.radial_factor
        expected = h^2 * radial_second_derivative(state[index], state[3 + index] / h,
                                                  mixed, 0, radial) + h_prime * state[3 + index]
        @test isapprox(derivative[3 + index], expected; rtol=1e-12)
    end
end
