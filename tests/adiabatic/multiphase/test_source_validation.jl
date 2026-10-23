@testset "Complete two-fundamental quartic table" begin
    grid = angular_grid(16, 9)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    phase1, phase2 = grid.azimuth_nodes, 2pi .* (0:10) ./ 11
    first_coefficient, second_coefficient = 0.2 + 0.1im, -0.07 + 0.04im
    radial = Dict((1, (1, 0)) => first_coefficient,
                  (2, (0, 1)) => second_coefficient)
    first_pattern = first_coefficient .* polar[(1, 1)]
    second_pattern = second_coefficient .* polar[(2, 2)]
    # Exact ordered-triple multiplicities from section 3.5, at each polar node.
    expected = Dict(
        (1, (1, 0)) => (abs2.(first_pattern)/2 + abs2.(second_pattern)) .* first_pattern,
        (2, (0, 1)) => (abs2.(second_pattern)/2 + abs2.(first_pattern)) .* second_pattern,
        (0, (2, -1)) => first_pattern.^2 .* conj.(second_pattern)/2,
        (3, (-1, 2)) => second_pattern.^2 .* conj.(first_pattern)/2,
        (3, (3, 0)) => first_pattern.^3/6,
        (6, (0, 3)) => second_pattern.^3/6,
        (3, (1, 1)) => zero.(first_pattern),
    )
    field = reconstruct_two_phase_field_at_radius(radial, [1, 2], polar, phase1, phase2)
    remainder = potential_remainder.(field, Ref(:quartic))
    radius, spin, alpha = 8.0, 0.7, 0.2
    source = nonlinear_two_phase_source(radial, radius, spin, alpha, keys(expected),
        [1, 2], polar, grid, phase1, phase2; potential=:quartic)
    for ((l, q), phase_coefficient) in expected
        recovered = [sum(remainder[i, j, k] * cis(q[1]*phase1[j] + q[2]*phase2[k])
            for j in eachindex(phase1), k in eachindex(phase2)) /
            (length(phase1)*length(phase2)) for i in eachindex(grid.polar_weights)]
        @test isapprox(recovered, phase_coefficient; atol=1e-13, rtol=1e-12)
        m = q[1] + 2q[2]
        integrated = alpha^2 * sum(grid.polar_weights .* polar[(l, abs(m))] .*
            (radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2) .* phase_coefficient)
        @test isapprox(source[(l, q)], integrated; atol=1e-13, rtol=1e-12)
    end
end

# A fixed, appreciably nonlinear two-cloud field. Keep input/output supports
# fixed while varying quadrature; omitted-source values do not solve profiles.
function cosine_source_grid_sample(npolar, nfirst, nsecond)
    grid = angular_grid(npolar, nfirst)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    phase2 = 2pi .* (0:nsecond-1) ./ nsecond
    radial = Dict((1, (1, 0)) => 2.3 + 0.7im, (2, (0, 1)) => 1.6 - 0.5im)
    outputs = [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2)),
               (3, (3, 0)), (6, (0, 3)), (5, (1, 0))]
    source = nonlinear_two_phase_source(radial, 8.0, 0.7, 0.2, outputs,
        [1, 2], polar, grid, grid.azimuth_nodes, phase2; potential=:cosine)
    return [source[pair] for pair in outputs]
end

@testset "Cosine polar and phase quadrature convergence" begin
    reference = cosine_source_grid_sample(64, 64, 64)
    for axis in 1:3
        errors = Float64[]
        for count in (8, 16, 32)
            counts = [64, 64, 64]
            counts[axis] = count
            source = cosine_source_grid_sample(counts...)
            push!(errors, maximum(abs.(source - reference)))
        end
        # Establish an actual convergence trend, rather than comparing two
        # grids that happen to alias to the same answer.
        @test errors[1] > 1e-5
        @test errors[2] < errors[1]/1000
        @test errors[3] < errors[2]
        @test errors[3] < 1e-10
        println("Cosine source axis ", axis, " errors at 8/16/32 nodes: ", errors)
    end
    cosine_source_grid_sample(32, 32, 32) # Warm this grid before measuring.
    measurement = @timed cosine_source_grid_sample(32, 32, 32)
    @test maximum(abs.(measurement.value - reference)) < 1e-10
    # Allocation count includes tables and temporary arrays; it is not peak RSS.
    println("Cosine 32^3 grid including setup: ", measurement.time,
        " seconds, ", measurement.bytes, " allocated bytes; one real field: ",
        32^3 * sizeof(Float64), " bytes")
end

@testset "Independent phase-origin shifts" begin
    grid = angular_grid(16, 32)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    phase1, phase2 = grid.azimuth_nodes, 2pi .* (0:35) ./ 36
    radial = Dict((1, (1, 0)) => 0.9 + 0.3im, (2, (0, 1)) => 0.4 - 0.2im)
    outputs = [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)),
               (3, (-1, 2)), (3, (3, 0)), (6, (0, 3))]
    field = reconstruct_two_phase_field_at_radius(radial, [1, 2], polar, phase1, phase2)
    for potential in (:quartic, :sextic, :cosine)
        source = nonlinear_two_phase_source(radial, 8.0, 0.7, 0.2, outputs,
            [1, 2], polar, grid, phase1, phase2; potential)
        average_potential = sum(potential_value.(field, Ref(potential)); dims=(2, 3)) /
            (length(phase1)*length(phase2))
        # Neither shift is an integer grid translation. Check each phase, then both.
        for shift in ((0.37, 0.0), (0.0, -0.53), (0.37, -0.53))
            shifted_radial = Dict((l, q) => value * cis(-sum(q .* shift))
                for ((l, q), value) in radial)
            shifted_source = nonlinear_two_phase_source(shifted_radial, 8.0, 0.7, 0.2,
                outputs, [1, 2], polar, grid, phase1, phase2; potential)
            @test isapprox([shifted_source[pair] for pair in outputs],
                [source[(l, q)] * cis(-sum(q .* shift)) for (l, q) in outputs];
                atol=1e-12, rtol=1e-12)
            shifted_field = reconstruct_two_phase_field_at_radius(
                shifted_radial, [1, 2], polar, phase1, phase2)
            shifted_average = sum(potential_value.(shifted_field, Ref(potential)); dims=(2, 3)) /
                (length(phase1)*length(phase2))
            # Invariance holds separately at every polar node, before polar integration.
            @test isapprox(shifted_average, average_potential; atol=1e-12, rtol=1e-12)
        end
    end
end
