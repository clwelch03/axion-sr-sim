@testset "Azimuthal average reproduces sum of |F_m|^2" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)

    # Arbitrary amplitudes on the primary's mode set, all three m blocks present.
    radial_values = Dict(
        (1, 1) => 0.3 + 0.2im,
        (3, 1) => -0.05 + 0.1im,
        (5, 1) => 0.01 - 0.02im,
        (3, 3) => 0.02 - 0.03im,
        (5, 3) => -0.004im,
        (5, 5) => 0.01im,
    )
    field = reconstruct_field(radial_values, polar, azimuthal, grid)

    # F_m(x) = Σ_l R_lm Y_lm(x) at each polar node.
    angular_profile = Dict(
        m => sum(amplitude .* polar[(l, mode_m)]
                 for ((l, mode_m), amplitude) in radial_values if mode_m == m)
        for m in (1, 3, 5)
    )
    expected = sum(abs2.(angular_profile[m]) for m in (1, 3, 5))

    # phi^2/2 has harmonics up to m = 10 < 31, so the average is exact.
    @test isapprox(azimuthal_average(field .^ 2 ./ 2, grid), expected;
                   atol=1e-14, rtol=1e-12)
end

@testset "Background coupling with unit deficit equals r^2 + a^2 cos^2 couplings" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(6, polar, grid)
    secondary_modes = [(2, 2), (4, 2), (6, 2)]
    radius, spin = 7.3, 0.99

    unit_deficit = ones(length(grid.cosine_theta_nodes))
    coupling_matrix = background_coupling_matrix(
        unit_deficit, radius, spin, secondary_modes, polar, grid,
    )

    expected = [
        radius^2 * (l == other_l) + spin^2 * couplings[(l, other_l, 2)]
        for (l, _) in secondary_modes, (other_l, _) in secondary_modes
    ]
    @test isapprox(coupling_matrix, expected; atol=1e-12, rtol=0)
    @test isapprox(coupling_matrix, transpose(coupling_matrix))
end

@testset "Fixed-time background term: direct part equals the period average" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6, grid.azimuth_nodes)
    radius, spin, alpha = 7.3, 0.99, 0.42

    # Order-one primary field, so the e^{4iφ} part of 1 - cos ϕ₁ is not small.
    primary_values = Dict(
        (1, 1) => 0.8 + 0.5im, (3, 1) => -0.1 + 0.2im, (5, 1) => 0.02im,
        (3, 3) => 0.15 - 0.1im, (5, 3) => 0.03 + 0.0im, (5, 5) => -0.02 + 0.01im,
    )
    secondary_modes = [(2, 2), (4, 2), (6, 2)]
    secondary_values = Dict((2, 2) => 0.3 - 0.1im, (4, 2) => 0.05 + 0.02im, (6, 2) => -0.01im)
    rotated_values = Dict(mode => im * value for (mode, value) in secondary_values)

    source = secondary_background_source(primary_values, secondary_values,
        radius, spin, alpha, secondary_modes, polar, azimuthal, grid)
    rotated_source = secondary_background_source(primary_values, rotated_values,
        radius, spin, alpha, secondary_modes, polar, azimuthal, grid)

    # Period-averaged term, computed directly.
    deficit = averaged_mass_deficit(primary_values, polar, azimuthal, grid)
    rho_squared = radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2
    secondary_profile = sum(value .* polar[mode] for (mode, value) in secondary_values)

    for mode in secondary_modes
        averaged = alpha^2 * sum(grid.polar_weights .* rho_squared .* deficit .*
                                 secondary_profile .* polar[mode])
        direct_part = (source[mode] - im * rotated_source[mode]) / 2
        @test isapprox(direct_part, averaged; atol=1e-12, rtol=1e-10)
    end
end

@testset "Secondary RHS with zero primary equals the linear m = 2 RHS" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(6, polar, grid)
    spin, alpha, frequency = 0.99, 0.42, 0.415 + 1e-7im
    secondary_modes = [(2, 2), (4, 2), (6, 2)]
    state = ComplexF64[0.3 - 0.1im, 0.05 + 0.02im, -0.01im, 0.02, -0.004im, 0.001 + 0.0im]

    zero_primary = Dict(mode => 0.0im for mode in [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)])

    secondary = similar(state)
    secondary_tortoise_rhs!(secondary, state, 7.3, zero_primary, spin, alpha, frequency,
                            secondary_modes, couplings, polar, azimuthal, grid)
    linear = similar(state)
    linear_tortoise_rhs!(linear, state, 7.3, spin, alpha, frequency,
                         secondary_modes, couplings; fundamental_m=2)

    @test isapprox(secondary, linear; atol=1e-14, rtol=0)
end

@testset "Secondary matching with zero primary equals linear m = 2 matching" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(6, polar, grid)
    spin, alpha, frequency = 0.99, 0.42, 0.415 + 1e-7im
    primary_modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
    secondary_modes = [(2, 2), (4, 2), (6, 2)]

    # A short synthetic primary domain with zero field, spacing 0.1.
    inner, match, outer, primary_spacing = -20.0, 5.0, 40.0, 0.1
    function zero_branch(start, stop)
        sample_count = round(Int, abs(stop - start) / primary_spacing) + 1
        return (; positions=collect(range(start, stop; length=sample_count)),
                  states=zeros(ComplexF64, 2length(primary_modes), sample_count))
    end
    zero_profiles = (; from_horizon=zero_branch(inner, match),
                       from_outer=zero_branch(outer, match))

    horizon_amplitudes = ComplexF64[1.0, 0.1im, -0.02]
    outer_amplitudes = ComplexF64[1.0, 0.05, 0.01im]

    secondary = secondary_matching_residual(
        frequency, horizon_amplitudes, outer_amplitudes, zero_profiles, primary_modes,
        spin, alpha, secondary_modes, couplings, polar, azimuthal, grid,
    )
    linear = matching_residual(
        frequency, horizon_amplitudes, outer_amplitudes, spin, alpha,
        [2, 4, 6], 2, couplings; step=2primary_spacing, inner, outer, match,
    )

    @test isapprox(secondary, linear; atol=1e-12, rtol=0)
end

@testset "Averaged background term: complex-linear and equal to the fixed-time direct part" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6, grid.azimuth_nodes)
    radius, spin, alpha = 7.3, 0.99, 0.42

    # Order-one primary field, so the fixed-time conjugate part is not small.
    primary_values = Dict(
        (1, 1) => 0.8 + 0.5im, (3, 1) => -0.1 + 0.2im, (5, 1) => 0.02im,
        (3, 3) => 0.15 - 0.1im, (5, 3) => 0.03 + 0.0im, (5, 5) => -0.02 + 0.01im,
    )
    secondary_modes = [(2, 2), (4, 2), (6, 2)]
    secondary_values = Dict((2, 2) => 0.3 - 0.1im, (4, 2) => 0.05 + 0.02im, (6, 2) => -0.01im)
    rotated_values = Dict(mode => im * value for (mode, value) in secondary_values)
    source(values; reading) = secondary_background_source(primary_values, values,
        radius, spin, alpha, secondary_modes, polar, azimuthal, grid; reading)

    averaged = source(secondary_values; reading=:averaged)
    averaged_rotated = source(rotated_values; reading=:averaged)
    fixed_time = source(secondary_values; reading=:fixed_time)
    fixed_time_rotated = source(rotated_values; reading=:fixed_time)

    for mode in secondary_modes
        # No conj(R): multiplying R by i multiplies the term by i.
        @test isapprox(averaged_rotated[mode], im * averaged[mode]; atol=1e-14, rtol=1e-12)
        # The fixed-time term's part proportional to R.
        direct_part = (fixed_time[mode] - im * fixed_time_rotated[mode]) / 2
        @test isapprox(averaged[mode], direct_part; atol=1e-12, rtol=1e-10)
    end
end

@testset "Uniform quadrature integrates cubics exactly, odd or even intervals" begin
    for sample_count in (8, 9)   # 7 and 8 intervals
        positions = range(0.0, 2.0; length=sample_count)
        @test isapprox(integrate_uniform(positions .^ 3, step(positions)), 4.0; atol=1e-13)
    end
end

