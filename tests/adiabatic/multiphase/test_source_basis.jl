@testset "Required source basis and direct quadrature agreement" begin
    grid = angular_grid(8, 9)
    polar_values = build_polar_values(3, grid.cosine_theta_nodes)
    phase1_nodes = grid.azimuth_nodes .+ 0.37
    phase2_nodes = 2pi .* (0:10) ./ 11 .- 0.53
    cloud_m = [1, 2]
    mode_labels = [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]
    source_basis = build_two_phase_source_basis(mode_labels, cloud_m,
        polar_values, grid, phase1_nodes, phase2_nodes)
    num_radial_modes = length(mode_labels)
    sample_coefficients = ComplexF64[0.2+0.1im, -0.07+0.04im, 0.01-0.02im, 0.03im]

    # The former direct-quadrature calculation remains a test reference only.
    # It reconstructs and projects with explicit phase sums, independently of
    # the precomputed matrix operations used by the solver.
    function direct_source(parameters; potential)
        radial_coefficients = complex.(parameters[1:num_radial_modes],
            parameters[num_radial_modes+1:2num_radial_modes])
        radius, spin, alpha = parameters[end-2:end]
        field = reconstruct_two_phase_field_at_radius(
            Dict(zip(mode_labels, radial_coefficients)), cloud_m, polar_values,
            phase1_nodes, phase2_nodes)
        geometry_factor = reshape(radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2, :, 1, 1)
        projected = two_phase_field_component_coefficients(
            geometry_factor .* potential_remainder.(field, Ref(potential)),
            mode_labels, cloud_m, polar_values, grid.polar_weights,
            phase1_nodes, phase2_nodes)
        return alpha^2 .* pack_complex([projected[label] for label in mode_labels])
    end

    for amplitude in (0.0, 1e-4, 1.0), potential in (:cosine, :chiral, :quartic, :sextic)
        parameters = vcat(pack_complex(amplitude .* sample_coefficients), [8.0, 0.7, 0.2])
        reference = values -> direct_source(values; potential)
        basis_source = values -> nonlinear_two_phase_source(values, source_basis; potential)
        @test isapprox(basis_source(parameters), reference(parameters); rtol=1e-11, atol=1e-24)
        reference_jacobian = ForwardDiff.jacobian(reference, parameters)
        @test isapprox(ForwardDiff.jacobian(basis_source, parameters), reference_jacobian;
            rtol=1e-10, atol=1e-20)
        local_source, local_jacobian = local_source_jacobian(parameters, source_basis; potential)
        @test isapprox(local_source, reference(parameters); rtol=1e-11, atol=1e-24)
        @test isapprox(local_jacobian, reference_jacobian; rtol=1e-10, atol=1e-20)
    end

    parameters = vcat(pack_complex(sample_coefficients), [8.0, 0.7, 0.2])
    @test_throws DimensionMismatch nonlinear_two_phase_source(parameters[1:end-1], source_basis)
    @test_throws DimensionMismatch nonlinear_two_phase_source(vcat(parameters, 0.0), source_basis)
end
