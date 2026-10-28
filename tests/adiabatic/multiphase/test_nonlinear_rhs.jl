# Compare the library RHS with the existing single-cloud implementation and
# exact driven-source formulas, including state ordering and tortoise scaling.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "Two-phase RHS single-cloud limit and ordering" begin
    grid = angular_grid(24, 64)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    phase1, phase2 = -grid.azimuth_nodes, 2pi .* (0:6) ./ 7
    cloud_m = [1, 2]
    frequencies = [0.405 + 1e-7im, 0.415 + 2e-7im]
    components = [fourier_component((1, 0), [1, 3], cloud_m),
                  fourier_component((3, 0), [3, 5], cloud_m)]
    modes = [(1, 1), (3, 1), (3, 3), (5, 3)]
    # Populate angular mixing and a higher harmonic, as well as nonzero slopes.
    state = ComplexF64[0.8+0.3im, -0.1+0.04im, 0.07-0.02im, -0.03+0.01im,
                       0.02+0.01im, -0.01+0.02im, 0.03-0.01im, -0.02-0.01im]
    original = copy(state)
    position, spin, alpha = 4.0, 0.7, 0.42
    mode_labels = [(l, component.q) for component in components for l in component.l_values]
    source_basis = build_two_phase_source_basis(mode_labels, cloud_m, polar, grid, phase1, phase2)
    for potential in (:quartic, :sextic, :cosine)
        derivative, reference = similar(state), similar(state)
        returned = nonlinear_two_phase_rhs!(derivative, state, position, spin, alpha,
            components, frequencies, cloud_m, couplings, source_basis;
            potential)
        nonlinear_tortoise_rhs!(reference, state, position, spin, alpha,
            frequencies[1], modes, couplings, polar, azimuthal, grid; potential)
        @test returned === derivative
        @test state == original
        @test derivative[1:4] == state[5:8]
        @test isapprox(derivative, reference; atol=1e-12, rtol=1e-12)

        # Reverse both component order and angular order inside each block.
        reordered = [fourier_component((3, 0), [5, 3], cloud_m),
                     fourier_component((1, 0), [3, 1], cloud_m)]
        permutation = [4, 3, 2, 1, 8, 7, 6, 5]
        reordered_basis = build_two_phase_source_basis(mode_labels[[4, 3, 2, 1]],
            cloud_m, polar, grid, phase1, phase2)
        reordered_derivative = similar(state)
        nonlinear_two_phase_rhs!(reordered_derivative, state[permutation], position,
            spin, alpha, reordered, frequencies, cloud_m, couplings, reordered_basis; potential)
        @test isapprox(reordered_derivative, derivative[permutation]; atol=1e-12, rtol=1e-12)
        @test_throws ArgumentError nonlinear_two_phase_rhs!(reordered_derivative,
            state[permutation], position, spin, alpha, reordered, frequencies,
            cloud_m, couplings, source_basis; potential)
    end
    @test_throws DimensionMismatch nonlinear_two_phase_rhs!(zeros(ComplexF64, 7),
        state, position, spin, alpha, components, frequencies, cloud_m, couplings,
        source_basis)
    @test_throws DimensionMismatch nonlinear_two_phase_rhs!(similar(state),
        state[1:7], position, spin, alpha, components, frequencies, cloud_m, couplings,
        source_basis)
end

@testset "Driven RHS sign and tortoise scaling" begin
    grid = angular_grid(16, 9)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar, grid)
    phase1, phase2 = grid.azimuth_nodes, 2pi .* (0:10) ./ 11
    cloud_m, frequencies = [1, 2], [0.405, 0.415]
    components = [fourier_component(q, [l], cloud_m)
        for (l, q) in [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]]
    mode_labels = [(l, component.q) for component in components for l in component.l_values]
    source_basis = build_two_phase_source_basis(mode_labels, cloud_m, polar, grid, phase1, phase2)
    first_coefficient, second_coefficient = 0.2+0.1im, -0.07+0.04im
    state = ComplexF64[first_coefficient, second_coefficient, 0, 0, 0, 0, 0, 0]
    position, spin, alpha = 4.0, 0.7, 0.42
    geometry = tortoise_geometry(position, spin)
    radius = geometry.radius
    source_scale = (radius^2 - 2radius + spin^2)/(radius^2 + spin^2)^2
    expected_first = -source_scale * alpha^2 * first_coefficient^2 *
        conj(second_coefficient) * sqrt(30)/20 * (radius^2 + spin^2/7)
    expected_second = -source_scale * alpha^2 * second_coefficient^2 *
        conj(first_coefficient) * sqrt(210)/42 * (radius^2 + spin^2/11)
    derivative = similar(state)
    nonlinear_two_phase_rhs!(derivative, state, position, spin, alpha, components,
        frequencies, cloud_m, couplings, source_basis; potential=:quartic)
    # R=Q=0 in the two driven blocks, so only their nonlinear source contributes.
    @test isapprox(derivative[7], expected_first; atol=1e-14, rtol=1e-12)
    @test isapprox(derivative[8], expected_second; atol=1e-14, rtol=1e-12)

    # With R=0 everywhere and nonzero Q, no source or V*R term can remain.
    slope_state = vcat(zeros(ComplexF64, 4), fill(0.03-0.02im, 4))
    nonlinear_two_phase_rhs!(derivative, slope_state, position, spin, alpha, components,
        frequencies, cloud_m, couplings, source_basis; potential=:quartic)
    @test derivative[1:4] == slope_state[5:8]
    @test isapprox(derivative[5:8], -2radius * source_scale .* slope_state[5:8];
        atol=1e-14, rtol=1e-12)
end
