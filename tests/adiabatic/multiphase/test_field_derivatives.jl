# Candidate checkpoint: load build_polar_theta_derivatives and
# reconstruct_two_phase_field_derivatives before this file.
# Register after the library accepts grid and builds the polar derivatives internally.
using ForwardDiff

@testset "Polar theta derivatives" begin
    grid = angular_grid(16, 9)
    polar = build_polar_values(4, grid.cosine_theta_nodes)
    derivatives = build_polar_theta_derivatives(polar, grid.cosine_theta_nodes)
    for (l, m) in keys(polar)
        reference = [ForwardDiff.derivative(theta -> polar_harmonic(l, m, cos(theta)),
            acos(x)) for x in grid.cosine_theta_nodes]
        @test isapprox(derivatives[(l, m)], reference; atol=1e-12, rtol=1e-12)
    end
    @test derivatives[(1, 1)] ≈ -sqrt(3/4) .* grid.cosine_theta_nodes
    @test all(iszero, derivatives[(0, 0)])
    @test_throws ArgumentError build_polar_theta_derivatives(polar, [-1.0, 1.0])
end

@testset "Frozen two-phase field and coordinate derivatives" begin
    grid = angular_grid(8, 9)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    phase1_nodes, phase2_nodes = grid.azimuth_nodes, 2pi .* (0:10) ./ 11
    position, spin = 8.0, 0.7
    frequencies, cloud_m = [0.4, 0.41], [1, 2]
    # Includes two distinct phase labels sharing l=m=1, m=0, negative m,
    # and two angular terms in one Fourier block.
    radial_values = Dict((1, (1, 0)) => 0.3+0.2im,
        (3, (1, 0)) => 0.04-0.03im, (2, (0, 1)) => -0.2+0.1im,
        (0, (2, -1)) => 0.1+0.06im, (1, (3, -2)) => -0.07+0.05im,
        (1, (3, -1)) => 0.09-0.08im)
    slopes = Dict(pair => (0.02-0.03im)*value for (pair, value) in radial_values)
    original_values, original_slopes = copy(radial_values), copy(slopes)
    reconstruct = (values, derivatives, carriers=frequencies) ->
        reconstruct_two_phase_field_derivatives(values, derivatives, position,
            spin, carriers, cloud_m, polar, grid, phase1_nodes, phase2_nodes)
    result = reconstruct(radial_values, slopes)
    geometry = tortoise_geometry(position, spin)
    @test size(result.field) == (8, 9, 11)
    @test result.field == reconstruct_two_phase_field_at_radius(radial_values,
        cloud_m, polar, phase1_nodes, phase2_nodes)

    # Differentiate an independently written real field in physical coordinates.
    # Each radial coefficient is locally affine in r, with prescribed dR/dr_*.
    maximum_errors = zeros(4)
    for i in (2, 6), j in (2, 5, 8), k in (3, 9)
        function explicit_field(coordinates)
            time, radius, theta, azimuth = coordinates
            sum(begin
                omega, m = sum(q .* frequencies), sum(q .* cloud_m)
                angle = q[1]*phase1_nodes[j] + q[2]*phase2_nodes[k] +
                    omega*time - m*azimuth
                value = coefficient + slopes[(l, q)] *
                    (radius-geometry.radius) / geometry.dr_drstar
                2polar_harmonic(l, abs(m), cos(theta)) *
                    (real(value)*cos(angle) + imag(value)*sin(angle))
            end for ((l, q), coefficient) in radial_values)
        end
        reference = ForwardDiff.gradient(explicit_field,
            [0.0, geometry.radius, acos(grid.cosine_theta_nodes[i]), 0.0])
        computed = [result.derivative_time[i,j,k], result.derivative_radius[i,j,k],
            result.derivative_theta[i,j,k], result.derivative_azimuth[i,j,k]]
        maximum_errors .= max.(maximum_errors, abs.(computed .- reference))
    end
    for error in maximum_errors
        @test error < 1e-12
    end
    @test radial_values == original_values && slopes == original_slopes
    zero_result = reconstruct(Dict(), Dict())
    @test all(array -> all(iszero, array), values(zero_result))
    axisymmetric = reconstruct(Dict((0, (2, -1)) => 0.1+0.06im),
        Dict((0, (2, -1)) => 0.0im))
    @test all(iszero, axisymmetric.derivative_azimuth)
    @test maximum(abs, axisymmetric.derivative_time) > 0.01
    @test_throws ArgumentError reconstruct(radial_values, slopes, complex.(frequencies))
    @test_throws ArgumentError reconstruct(radial_values, slopes, [NaN, 0.41])
    @test_throws DimensionMismatch reconstruct(radial_values, slopes, [0.4])
    @test_throws ArgumentError reconstruct(radial_values, Dict())
    println("Maximum t/r/theta/azimuth derivative errors: ", maximum_errors)
end
