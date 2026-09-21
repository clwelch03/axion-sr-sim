using LinearAlgebra: norm

"""
    differentiated_series(coefficients, argument) -> (value, first, second)

Evaluate a vector polynomial and its first two derivatives term by term.

Store coefficients in ascending powers. Differentiate with respect to `argument`
to provide an independent reference for boundary ODE residual tests.
"""
function differentiated_series(coefficients, argument)
    value, first, second = zero(coefficients[1]), zero(coefficients[1]), zero(coefficients[1])
    for (index, coefficient) in enumerate(coefficients)
        power = index - 1
        value += coefficient * argument^power
        power >= 1 && (first += power * coefficient * argument^(power - 1))
        power >= 2 && (second += power * (power - 1) * coefficient * argument^(power - 2))
    end
    return value, first, second
end

@testset "Boundary sheets and asymptotic signs" begin
    for growth in (-1e-7, 0.0, 1e-7)
        bound = outer_wave_number(0.4088 + growth * im, 0.42; kind=:bound)
        outgoing = outer_wave_number(1.2 + growth * im, 0.42; kind=:outgoing)
        @test imag(bound) > 0
        @test real(outgoing) > 0
        @test bound^2 ≈ (0.4088 + growth * im)^2 - 0.42^2
        @test outgoing^2 ≈ (1.2 + growth * im)^2 - 0.42^2
    end
    @test_throws DomainError outer_wave_number(0.42, 0.42; kind=:outgoing)
    @test_throws DomainError outer_wave_number(1.2, 0.42; kind=:bound)
    @test_throws ArgumentError outer_wave_number(0.4, 0.42; kind=:unknown)
end

@testset "Boundary data and derivatives" begin
    grid = angular_grid(16, 31)
    couplings = build_cosine_squared_couplings(5, build_polar_values(5, grid.cosine_theta_nodes), grid)
    l_values, amplitudes = [1, 3, 5], ComplexF64[1, 0.1im, -0.05]
    spin, alpha, frequency = 0.99, 0.42, 0.4088 + 1.5e-7im
    # A natural amplitude input [1,0,0] must not force integer series coefficients.
    for boundary in (horizon_series, outer_series)
        integer_input = boundary(spin, alpha, frequency, l_values, 1, couplings, [1, 0, 0])
        complex_input = boundary(spin, alpha, frequency, l_values, 1, couplings, ComplexF64[1, 0, 0])
        @test integer_input.coefficients ≈ complex_input.coefficients
    end
    inner = horizon_boundary_data(-400.0, spin, alpha, frequency, l_values, 1, couplings, amplitudes)
    horizon_frequency = spin / (2kerr_horizons(spin).r_plus)
    @test inner.tortoise_slopes ≈ -im * (frequency - horizon_frequency) * inner.values
    outer = outer_boundary_data(100.0, spin, alpha, frequency, l_values, 1, couplings, amplitudes; order=0)
    geometry = tortoise_geometry(100.0, spin)
    decay = sqrt(complex(alpha^2 - frequency^2))
    power = -1 - (alpha^2 - 2frequency^2) / decay
    @test outer.values ≈ exp(-decay * geometry.radius) * geometry.radius^power * amplitudes
    @test outer.tortoise_slopes ≈ geometry.dr_drstar * (-decay + power / geometry.radius) * outer.values
    # Derivatives must be with respect to r_*, not r.
    for (boundary, coordinate) in ((horizon_boundary_data, -40.0), (outer_boundary_data, 100.0))
        step = 1e-4
        data = boundary(coordinate, spin, alpha, frequency, l_values, 1, couplings, amplitudes)
        plus = boundary(coordinate + step, spin, alpha, frequency, l_values, 1, couplings, amplitudes)
        minus = boundary(coordinate - step, spin, alpha, frequency, l_values, 1, couplings, amplitudes)
        @test isapprox((plus.values - minus.values) / (2step), data.tortoise_slopes; rtol=2e-7)
    end
end

@testset "Asymptotic series satisfy the coupled radial equation" begin
    grid = angular_grid(16, 31)
    couplings = build_cosine_squared_couplings(5, build_polar_values(5, grid.cosine_theta_nodes), grid)
    l_values, amplitudes = [1, 3, 5], ComplexF64[1, 0.1im, -0.05]
    spin, alpha = 0.99, 0.42
    for (frequency, kind) in ((0.4088 + 1.5e-7im, :bound), (1.2 + 1e-7im, :outgoing))
        outer_residuals, inner_residuals = Float64[], Float64[]
        for order in (0, 3, 7)
            series = outer_series(spin, alpha, frequency, l_values, 1, couplings, amplitudes; kind, order)
            radius = 200.0
            value, first, second = differentiated_series(series.coefficients, 1 / radius)
            logarithmic_slope = im * series.wave_number + series.exponent / radius
            radial_first = logarithmic_slope * value - first / radius^2
            radial_second = (logarithmic_slope^2 - series.exponent / radius^2) * value -
                2logarithmic_slope * first / radius^2 + 2first / radius^3 + second / radius^4
            residual = similar(value)
            for (index, l) in enumerate(l_values)
                coefficients = radial_operator_coefficients(radius, spin, alpha, frequency, l, 1)
                mixed = sum(couplings[(l, other_l, 1)] * value[j] for (j, other_l) in enumerate(l_values))
                residual[index] = coefficients.delta * radial_second[index] + coefficients.delta_prime * radial_first[index] +
                    coefficients.local_coefficient * value[index] + coefficients.coupling_scale * mixed
            end
            push!(outer_residuals, norm(residual) / (radius^2 * norm(value)))

            series = horizon_series(spin, alpha, frequency, l_values, 1, couplings, amplitudes; order)
            distance = 0.001
            radius = kerr_horizons(spin).r_plus + distance
            gap = kerr_horizons(spin).gap
            value, first, second = differentiated_series(series.coefficients, distance)
            p = series.exponent
            scaled_first = p * value + distance * first
            scaled_second = p * (p - 1) * value + 2p * distance * first + distance^2 * second
            for (index, l) in enumerate(l_values)
                coefficients = radial_operator_coefficients(radius, spin, alpha, frequency, l, 1)
                mixed = sum(couplings[(l, other_l, 1)] * value[j] for (j, other_l) in enumerate(l_values))
                residual[index] = (distance + gap)^2 * scaled_second[index] +
                    (distance + gap) * (2distance + gap) * scaled_first[index] +
                    coefficients.delta * (coefficients.local_coefficient * value[index] + coefficients.coupling_scale * mixed)
            end
            push!(inner_residuals, norm(residual) / norm(value))
        end
        @test outer_residuals[2] < outer_residuals[1]
        @test outer_residuals[3] < outer_residuals[2]
        @test outer_residuals[3] < 1e-10
        @test inner_residuals[2] < inner_residuals[1]
        @test inner_residuals[3] < 1e-10
        println("$kind boundary residuals, orders 0/3/7: outer=$outer_residuals; horizon=$inner_residuals")
    end
end
