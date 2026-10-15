# Candidate checkpoint: load kerr_metric and scalar_stress_tensor first.
# Pointwise contractions below are reference test plumbing, not a library API.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))
using LinearAlgebra: Diagonal, issymmetric

@testset "Scalar stress tensor in flat spacetime" begin
    metric = Matrix(Diagonal([-1.0, 1.0, 1.0, 1.0]))
    field, alpha = 0.8, 0.42
    gradient = [0.7, -0.2, 0.3, -0.4]
    for potential in (:cosine, :quartic, :sextic)
        tensor = scalar_stress_tensor(field, gradient, metric, metric, alpha; potential)
        expected_energy = sum(abs2, gradient)/2 + alpha^2*potential_value(field, potential)
        @test tensor[1, 1] ≈ expected_energy
        @test tensor[1, 2:4] ≈ gradient[1] .* gradient[2:4]
        @test issymmetric(tensor)
    end
    potential_only = scalar_stress_tensor(field, zeros(4), metric, metric, alpha)
    potential_energy = alpha^2 * (1-cos(field))
    @test potential_only ≈ Diagonal([potential_energy, -potential_energy,
        -potential_energy, -potential_energy])
    # A massless right-moving plane wave has equal energy density and +r energy flux.
    null_tensor = scalar_stress_tensor(0.0, [-2.0, 2.0, 0.0, 0.0], metric, metric, 0.0)
    @test null_tensor[1, 1] ≈ 4.0
    @test -(metric*null_tensor)[2, 1] ≈ 4.0
    @test_throws DimensionMismatch scalar_stress_tensor(field, gradient[1:3], metric, metric, alpha)
end

function stress_tensor_quadrature_for_test(samples, position, spin, alpha,
        grid; potential=:cosine)
    geometry = tortoise_geometry(position, spin)
    phase_count = size(samples.field, 2) * size(samples.field, 3)
    energy = angular_momentum = energy_flux = angular_momentum_flux = 0.0
    for i in eachindex(grid.cosine_theta_nodes)
        x = grid.cosine_theta_nodes[i]
        metric = kerr_metric(position, spin, x)
        angular_weight = 2pi * grid.polar_weights[i] *
            (geometry.radius^2 + spin^2*x^2) / phase_count
        for j in axes(samples.field, 2), k in axes(samples.field, 3)
            gradient = [samples.derivative_time[i,j,k], samples.derivative_radius[i,j,k],
                samples.derivative_theta[i,j,k], samples.derivative_azimuth[i,j,k]]
            tensor = scalar_stress_tensor(samples.field[i,j,k], gradient,
                metric.metric, metric.inverse_metric, alpha; potential)
            mixed = metric.inverse_metric * tensor
            energy -= angular_weight * geometry.dr_drstar * mixed[1, 1]
            angular_momentum += angular_weight * geometry.dr_drstar * mixed[1, 4]
            energy_flux -= angular_weight * mixed[2, 1]
            angular_momentum_flux += angular_weight * mixed[2, 4]
        end
    end
    return (; energy, angular_momentum, energy_flux, angular_momentum_flux)
end

@testset "Kerr tensor energy agrees with single-cloud energy" begin
    grid = angular_grid(24, 32)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(3, grid.azimuth_nodes)
    modes = [(1, 1), (3, 1), (3, 3)]
    state = ComplexF64[0.7+0.2im, 0.1-0.03im, -0.04+0.02im,
        -0.03+0.01im, 0.004+0.002im, 0.002-0.001im]
    radial = Dict((l, (m, 0)) => state[index] for (index, (l, m)) in enumerate(modes))
    slopes = Dict((l, (m, 0)) => state[3+index] for (index, (l, m)) in enumerate(modes))
    # The second phase is inactive here. Both implementations use real carriers.
    position, alpha, frequency = 8.0, 0.42, 0.4
    errors = Float64[]
    for spin in (0.0, 0.7, 0.99), potential in (:cosine, :quartic, :sextic)
        samples = reconstruct_two_phase_field_derivatives(radial, slopes, position,
            spin, [frequency, 0.41], [1, 2], polar, grid, -grid.azimuth_nodes, [0.0])
        observables = stress_tensor_quadrature_for_test(samples, position, spin, alpha,
            grid; potential)
        reference = energy_per_tortoise_radius(state, position, frequency, spin, alpha,
            modes, polar, azimuthal, grid; potential)
        push!(errors, abs(observables.energy-reference)/abs(reference))
        @test isapprox(observables.energy, reference; rtol=2e-12)
    end
    println("Maximum single-cloud tensor energy relative error: ", maximum(errors))
end

@testset "Joint tensor angular momentum and channel flux signs" begin
    grid = angular_grid(16, 17)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    radial = Dict((1, (1, 0)) => 0.3+0.2im, (2, (0, 1)) => -0.2+0.1im,
        (0, (2, -1)) => 0.05-0.03im, (1, (3, -2)) => 0.03+0.01im)
    # Local outgoing-wave slopes: each channel has positive outward energy current.
    slopes = Dict(pair => 0.07im*value for (pair, value) in radial)
    frequencies, cloud_m = [0.4, 0.41], [1, 2]
    position, alpha = 8.0, 0.42
    for spin in (0.0, 0.7)
        samples = reconstruct_two_phase_field_derivatives(radial, slopes, position,
            spin, frequencies, cloud_m, polar, grid, grid.azimuth_nodes, grid.azimuth_nodes)
        result = stress_tensor_quadrature_for_test(samples, position, spin, alpha, grid)
        geometry = tortoise_geometry(position, spin)
        expected_energy_flux = 4pi * geometry.radial_factor * sum(
            sum(q .* frequencies) * imag(conj(value)*slopes[(l, q)])
            for ((l, q), value) in radial)
        expected_angular_flux = 4pi * geometry.radial_factor * sum(
            sum(q .* cloud_m) * imag(conj(value)*slopes[(l, q)])
            for ((l, q), value) in radial)
        @test result.energy_flux > 0
        @test isapprox(result.energy_flux, expected_energy_flux; rtol=2e-12)
        @test isapprox(result.angular_momentum_flux, expected_angular_flux; rtol=2e-12)
        if spin == 0.0
            expected_momentum = 4pi * geometry.radius^2 * sum(
                sum(q .* frequencies) * sum(q .* cloud_m) * abs2(value)
                for ((l, q), value) in radial)
            @test isapprox(result.angular_momentum, expected_momentum; rtol=2e-12)
        end
    end
end
