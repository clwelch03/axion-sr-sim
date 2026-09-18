@testset "Angular reconstruction and projection" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)

    amplitudes = Dict(
        (1, 1) => 0.3 + 0.2im,
        (2, 1) => -0.1 + 0.4im,  # Include an even-l mode.
        (3, 3) => 0.02 - 0.03im,
    )

    field = reconstruct_field(amplitudes, polar, azimuthal, grid)
    modes = [(l, m) for l in 1:5 for m in 1:l]
    recovered = project_field(field, modes, polar, azimuthal, grid)

    for mode in modes
        # Absent input modes should project to zero.
        expected = get(amplitudes, mode, 0.0im)
        @test isapprox(recovered[mode], expected; atol=1e-12, rtol=0)
    end

    # Independent phase check: R=i gives phi=-2Y*sin(varphi).
    imaginary_mode = reconstruct_field(
        Dict((1, 1) => 1im), polar, azimuthal, grid)
    expected_field = -2 .* polar[(1, 1)] * transpose(sin.(grid.azimuth_nodes))
    @test isapprox(imaginary_mode, expected_field; atol=1e-12, rtol=0)
end

@testset "Sine remainder near zero and across the evaluation switch" begin
    # BigFloat subtraction supplies an independent reference, including where
    # ordinary Float64 subtraction would erase the tiny cubic remainder.
    for field in (1e-12, 1e-6, 0.01, prevfloat(0.1), 0.1, nextfloat(0.1), 0.7, 3.0)
        reference = setprecision(BigFloat, 256) do
            exact_field = BigFloat(field)
            Float64(exact_field - sin(exact_field))
        end
        @test isapprox(sine_remainder(field), reference; rtol=2e-13, atol=0)
        @test sine_remainder(-field) == -sine_remainder(field)
    end
    @test sine_remainder(0.0) == 0
    @test isapprox(sine_remainder(1e-12) / (1e-12)^3, 1/6; rtol=1e-15)
end

@testset "Analytic cubic source and generated harmonics" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
    modes = [(l, m) for l in 1:5 for m in 1:l]
    radius, alpha, amplitude = 4.0, 0.42, 1e-4 * (1 + 0.4im)

    for spin in (0.0, 0.99)
        source = nonlinear_source(Dict((1, 1) => amplitude), radius, spin, alpha,
                                  polar, azimuthal, grid; modes)
        # phi^3/6 has m=1 coefficient |A|^2*A*Y11^3/2 and m=3
        # coefficient A^3*Y11^3/6. The polar integrals can be done exactly:
        # int Y11^4 dx=3/5, int x^2*Y11^4 dx=3/35;
        # int Y33*Y11^3 dx=sqrt(27/70), with x^2-weighted ratio 1/9.
        expected_11 = alpha^2 * abs2(amplitude) * amplitude * (3radius^2/10 + 3spin^2/70)
        expected_33 = alpha^2 * amplitude^3 * sqrt(3/280) * (radius^2 + spin^2/9)
        # O(|A|^2) relative corrections come from the full sine's fifth-order term.
        @test isapprox(source[(1, 1)], expected_11; rtol=1e-8, atol=0)
        @test isapprox(source[(3, 3)], expected_33; rtol=1e-8, atol=0)

        halved = nonlinear_source(Dict((1, 1) => amplitude/2), radius, spin, alpha,
                                  polar, azimuthal, grid; modes)
        for mode in ((1, 1), (3, 3))
            @test isapprox(source[mode], 8halved[mode]; rtol=1e-8, atol=0)
        end
        # A pure Y11 field is even in x and changes sign under varphi -> varphi+pi.
        # These output modes are therefore forbidden by polar or azimuthal parity.
        for mode in ((2, 1), (4, 1), (4, 3), (2, 2), (3, 2), (4, 4))
            @test abs(source[mode]) < 1e-12 * abs(source[(1, 1)])
        end
    end
    zero_source = nonlinear_source(Dict((1, 1) => 0.0im), radius, 0.99, alpha,
                                   polar, azimuthal, grid; modes)
    @test all(iszero, values(zero_source))
end

@testset "Full sine source against analytic angular moments" begin
    grid = angular_grid(32, 63)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
    radius, spin, alpha = 4.0, 0.99, 0.42
    for amplitude in (0.8, 1.5)
        # For real A and a pure Y11 mode, integrate the sine series analytically.
        # At order A^(2k+1), the azimuthal coefficient is 1/(k!*(k+1)!).
        # The polar moment is int Y11^(2p) dx, p=k+1. Its x^2-weighted
        # counterpart is smaller by 2p+3. No numerical quadrature enters here.
        reference = setprecision(BigFloat, 256) do
            polar_moment = BigFloat(1)  # int Y11^2 dx
            total = BigFloat(0)
            for order in 1:30
                power = order + 1
                polar_moment *= (big(3)/4) * (2power) / (2power + 1)
                coefficient = (-1)^(order + 1) * BigFloat(amplitude)^(2order + 1) /
                              (factorial(big(order)) * factorial(big(order + 1)))
                total += coefficient * polar_moment *
                         (BigFloat(radius)^2 + BigFloat(spin)^2 / (2power + 3))
            end
            Float64(BigFloat(alpha)^2 * total)
        end
        source = nonlinear_source(Dict((1, 1) => amplitude), radius, spin, alpha,
                                  polar, azimuthal, grid)
        @test isapprox(source[(1, 1)], reference; rtol=1e-12, atol=0)
    end
end

@testset "Multimode nonlinear projection: resolution and phase" begin
    modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
    amplitudes = Dict(zip(modes, ComplexF64[0.8+0.2im, -0.1im, 0.03,
                                           0.12-0.04im, 0.02im, -0.025]))
    function source_at_resolution(n_polar, n_azimuth, radial_values)
        grid = angular_grid(n_polar, n_azimuth)
        source = nonlinear_source(radial_values, 4.0, 0.99, 0.42,
            build_polar_values(5, grid.cosine_theta_nodes),
            build_azimuthal_values(5, grid.azimuth_nodes), grid; modes)
        return [source[mode] for mode in modes]
    end

    reference = source_at_resolution(64, 127, amplitudes)
    # Refine polar and azimuthal resolution independently, then together.
    for resolution in ((32, 127), (64, 63), (32, 63))
        source = source_at_resolution(resolution..., amplitudes)
        @test norm(source - reference) / norm(reference) < 1e-10
    end
    coarse = source_at_resolution(16, 31, amplitudes)
    println("Multimode nonlinear source: 16x31 vs 64x127 relative difference = ",
            norm(coarse - reference) / norm(reference))

    # An azimuthal shift rotates each amplitude and projected source by exp(i*m*chi).
    # The fine grid makes quadrature aliasing negligible in this symmetry check.
    phase = 0.37
    rotated = Dict((l, m) => amplitude * cis(m * phase)
                   for ((l, m), amplitude) in amplitudes)
    shifted_source = source_at_resolution(64, 127, rotated)
    expected = reference .* [cis(m * phase) for (l, m) in modes]
    @test norm(shifted_source - expected) / norm(reference) < 1e-12
end
