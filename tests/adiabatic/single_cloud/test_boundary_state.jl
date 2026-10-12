@testset "Boundary state harmonic frequencies and ordering" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(6, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(6, polar_values, grid)
    spin, alpha, frequency = 0.99, 0.42, 0.4088 + 1.5e-7im

    # Interleave m blocks so grouping must preserve both amplitude and state order.
    for (fundamental_m, modes) in (
        (1, [(3, 3), (5, 1), (5, 5), (1, 1), (5, 3), (3, 1)]),
        (2, [(6, 6), (4, 2), (2, 2)]),
    )
        count = length(modes)
        amplitudes = ComplexF64[0.1index + 0.03im * index^2 for index in 1:count]
        original_amplitudes = copy(amplitudes)
        outer_geometry = tortoise_geometry(100.0, spin)
        inner_geometry = tortoise_geometry(-40.0, spin)
        outer_state = boundary_state(:outer, 100.0, spin, alpha, frequency,
            modes, couplings, amplitudes; fundamental_m, outer_order=0)
        inner_state = boundary_state(:horizon, -40.0, spin, alpha, frequency,
            modes, couplings, amplitudes; fundamental_m, horizon_order=0)
        expected_outer = zeros(ComplexF64, 2count)
        expected_inner = similar(expected_outer)

        # Independent order-zero carriers: decaying for n=1, outgoing for n>1.
        for (index, (_, m)) in enumerate(modes)
            harmonic_frequency = (m ÷ fundamental_m) * frequency
            wave_number = m == fundamental_m ?
                im * sqrt(alpha^2 - harmonic_frequency^2) :
                sqrt(harmonic_frequency^2 - alpha^2)
            outer_power = -1 + im * (2harmonic_frequency^2 - alpha^2) / wave_number
            radius = outer_geometry.radius
            expected_outer[index] = amplitudes[index] *
                exp(im * wave_number * radius + outer_power * log(radius))
            expected_outer[count + index] = outer_geometry.dr_drstar *
                (im * wave_number + outer_power / radius) * expected_outer[index]

            horizon_power = -im * (2inner_geometry.r_plus * harmonic_frequency - spin * m) /
                inner_geometry.gap
            expected_inner[index] = amplitudes[index] *
                exp(horizon_power * inner_geometry.log_horizon_distance)
            expected_inner[count + index] = horizon_power *
                (inner_geometry.horizon_distance + inner_geometry.gap) /
                inner_geometry.radial_factor * expected_inner[index]
        end
        @test outer_state ≈ expected_outer
        @test inner_state ≈ expected_inner
        @test amplitudes == original_amplitudes

        # Include the coupled series corrections and verify permutation covariance.
        permutation = reverse(1:count)
        for (side, position) in ((:horizon, -40.0), (:outer, 100.0))
            state = boundary_state(side, position, spin, alpha, frequency,
                modes, couplings, amplitudes; fundamental_m)
            reordered = boundary_state(side, position, spin, alpha, frequency,
                modes[permutation], couplings, amplitudes[permutation]; fundamental_m)
            @test reordered ≈ state[vcat(permutation, count .+ permutation)]
        end
    end

    # With one m block, assembly must reproduce the existing boundary helpers.
    modes, amplitudes = [(1, 1), (3, 1), (5, 1)], ComplexF64[1, 0.1im, -0.05]
    for (side, position, boundary) in ((:horizon, -40.0, horizon_boundary_data),
                                      (:outer, 100.0, outer_boundary_data))
        reference = boundary(position, spin, alpha, frequency, [1, 3, 5], 1,
            couplings, amplitudes; order=2)
        state = boundary_state(side, position, spin, alpha, frequency,
            modes, couplings, amplitudes; horizon_order=2, outer_order=2)
        @test state ≈ vcat(reference.values, reference.tortoise_slopes)
    end

    @test_throws ArgumentError boundary_state(:unknown, 100.0, spin, alpha,
        frequency, modes, couplings, amplitudes)
    @test_throws DimensionMismatch boundary_state(:outer, 100.0, spin, alpha,
        frequency, modes, couplings, [1])
    for fundamental_m in (0, -1, 1.5, 2)
        @test_throws ArgumentError boundary_state(:outer, 100.0, spin, alpha,
            frequency, modes, couplings, amplitudes; fundamental_m)
    end
    @test_throws ArgumentError boundary_state(:horizon, -40.0, spin, alpha,
        frequency, [(0, 0)], couplings, [1])
end
