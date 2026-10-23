# Test library reconstruction and projection against known coefficients, keeping
# the explicit quadrature as a reference for the benchmark in plan section 3.3.
@testset "Independent phases separate equal angular patterns" begin
    grid = angular_grid(16, 9)
    polar_values = polar_harmonic.(3, 3, grid.cosine_theta_nodes)
    # Reuse the uniform periodic nodes for each independent phase coordinate.
    phase_nodes = grid.azimuth_nodes
    phase_count = length(phase_nodes)
    first_coefficient, second_coefficient = 0.2 + 0.1im, -0.07 + 0.04im

    # q=(3,0) and q=(-1,2) both have m=3 when cloud_m=(1,2).
    # Array axes are polar position, first phase, and second phase.
    radial_values = Dict((3, (3, 0)) => first_coefficient,
                         (3, (-1, 2)) => second_coefficient)
    polar_table = build_polar_values(3, grid.cosine_theta_nodes)
    field = reconstruct_two_phase_field_at_radius(radial_values, [1, 2],
        polar_table, phase_nodes, phase_nodes)

    recovered = ComplexF64[]
    for (q, expected) in (((3, 0), first_coefficient), ((-1, 2), second_coefficient))
        coefficient = sum(
            grid.polar_weights[i] * polar_values[i] *
            cis(q[1] * phase_nodes[j] + q[2] * phase_nodes[k]) * field[i, j, k]
            for i in eachindex(polar_values), j in eachindex(phase_nodes),
                k in eachindex(phase_nodes)
        ) / phase_count^2
        push!(recovered, coefficient)
        @test isapprox(coefficient, expected; atol=1e-13, rtol=0)
    end

    # Request both populated patterns, an absent pattern, and a conjugate.
    projected = two_phase_field_component_coefficients(field,
        [(3, (3, 0)), (3, (-1, 2)), (1, (1, 0)), (3, (-3, 0))],
        [1, 2], polar_table, grid.polar_weights, phase_nodes, phase_nodes)
    @test isapprox(projected[(3, (3, 0))], first_coefficient; atol=1e-13, rtol=0)
    @test isapprox(projected[(3, (-1, 2))], second_coefficient; atol=1e-13, rtol=0)
    @test isapprox(projected[(1, (1, 0))], 0; atol=1e-13, rtol=0)
    @test isapprox(projected[(3, (-3, 0))], conj(first_coefficient); atol=1e-13, rtol=0)

    # At zero temporal phases, psi1=-azimuth and psi2=-2azimuth.
    # Sample that slice of the SAME field, wrapping indices periodically.
    # An m=3 azimuthal projection sees the sum of the two coefficients.
    azimuthal_coefficient = sum(
        grid.polar_weights[i] * polar_values[i] * cis(-3phase_nodes[j + 1]) *
        field[i, mod(-j, phase_count) + 1, mod(-2j, phase_count) + 1]
        for i in eachindex(polar_values), j in 0:phase_count-1
    ) / phase_count
    @test isapprox(azimuthal_coefficient, first_coefficient + second_coefficient;
        atol=1e-13, rtol=0)

    println("Two-phase coefficient maximum error: ",
        maximum(abs.(recovered .- [first_coefficient, second_coefficient])))
end

@testset "Reconstruction shape, phase sign, and zero field" begin
    grid = angular_grid(16, 9)
    polar_table = build_polar_values(3, grid.cosine_theta_nodes)
    phase1_nodes = grid.azimuth_nodes
    phase2_nodes = 2pi .* (0:10) ./ 11

    field = reconstruct_two_phase_field_at_radius(Dict((1, (1, 0)) => 1im),
        [1, 2], polar_table, phase1_nodes, phase2_nodes)
    # 2Re[i exp(-i psi1)] = 2sin(psi1), independently of psi2.
    expected = [2polar*sin(phase1) for polar in polar_table[(1, 1)],
        phase1 in phase1_nodes, phase2 in phase2_nodes]
    @test size(field) == (16, 9, 11)
    @test field ≈ expected
    # Project the analytic field directly: test the sign and normalization on
    # unequal phase grids independently of the reconstruction implementation.
    projected = two_phase_field_component_coefficients(expected, [(1, (1, 0))], [1, 2],
        polar_table, grid.polar_weights, phase1_nodes, phase2_nodes)
    @test isapprox(projected[(1, (1, 0))], 1im; atol=1e-13, rtol=0)
    @test all(iszero, reconstruct_two_phase_field_at_radius(Dict(), [1, 2],
        polar_table, phase1_nodes, phase2_nodes))
    @test_throws DimensionMismatch reconstruct_two_phase_field_at_radius(
        Dict(), [1], polar_table, phase1_nodes, phase2_nodes)
    @test_throws DimensionMismatch reconstruct_two_phase_field_at_radius(
        Dict((1, (1,)) => 1), [1, 2], polar_table, phase1_nodes, phase2_nodes)
end
