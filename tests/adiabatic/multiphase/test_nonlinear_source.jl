# Load the single-cloud implementation only as a test reference.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

# Exact integrated quartic benchmarks from implementation-plan section 3.5.
# These expected values use analytic angular integrals, not a second projector.
@testset "Two-phase quartic nonlinear source" begin
    grid = angular_grid(16, 9)
    polar_table = build_polar_values(3, grid.cosine_theta_nodes)
    phase1_nodes = grid.azimuth_nodes
    phase2_nodes = 2pi .* (0:10) ./ 11
    first_coefficient, second_coefficient = 0.2 + 0.1im, -0.07 + 0.04im
    radial_values = Dict((1, (1, 0)) => first_coefficient,
                         (2, (0, 1)) => second_coefficient)
    output_pairs = [(0, (2, -1)), (3, (-1, 2)), (1, (2, -1)), (3, (1, 1))]
    radius, spin, alpha = 8.0, 0.7, 0.2

    mode_labels = union(collect(keys(radial_values)), output_pairs)
    source_basis = build_two_phase_source_basis(mode_labels, [1, 2], polar_table,
        grid, phase1_nodes, phase2_nodes)
    coefficients = [get(radial_values, label, 0.0im) for label in mode_labels]
    parameters = vcat(pack_complex(coefficients), [radius, spin, alpha])
    source = Dict(zip(mode_labels, unpack_complex(
        nonlinear_two_phase_source(parameters, source_basis; potential=:quartic))))

    # The cubic cross terms drive components absent from radial_values.
    expected_first = alpha^2 * first_coefficient^2 * conj(second_coefficient) *
        sqrt(30)/20 * (radius^2 + spin^2/7)
    expected_second = alpha^2 * second_coefficient^2 * conj(first_coefficient) *
        sqrt(210)/42 * (radius^2 + spin^2/11)
    @test isapprox(source[(0, (2, -1))], expected_first; atol=1e-14, rtol=1e-12)
    @test isapprox(source[(3, (-1, 2))], expected_second; atol=1e-14, rtol=1e-12)

    # Y_10 is odd in x, but the source at q=(2,-1) is even.
    @test isapprox(source[(1, (2, -1))], 0; atol=1e-14, rtol=0)
    # No sum of three signed input labels gives q=(1,1).
    @test isapprox(source[(3, (1, 1))], 0; atol=1e-14, rtol=0)

    empty_parameters = vcat(zeros(2length(mode_labels)), [radius, spin, alpha])
    empty_source = nonlinear_two_phase_source(empty_parameters, source_basis;
        potential=:quartic)
    @test all(iszero, empty_source)
end

@testset "Nonlinear source in the single-cloud limit" begin
    grid = angular_grid(24, 64)
    polar_table = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal_table = build_azimuthal_values(3, grid.azimuth_nodes)
    # At time zero for m_1=1, psi_1=-phi: exp(-i*j*psi_1)=exp(+i*j*phi).
    # Thus the same complex coefficient represents the same field in both APIs.
    phase1_nodes = -grid.azimuth_nodes
    phase2_nodes = 2pi .* (0:10) ./ 11
    radial_coefficient = 0.8 + 0.3im
    radial_values = Dict((1, (1, 0)) => radial_coefficient)
    single_cloud_values = Dict((1, 1) => radial_coefficient)
    output_pairs = [(1, (1, 0)), (3, (3, 0)), (3, (-1, 2))]
    radius, spin, alpha = 8.0, 0.7, 0.2
    source_basis = build_two_phase_source_basis(output_pairs, [1, 2], polar_table,
        grid, phase1_nodes, phase2_nodes)
    coefficients = [get(radial_values, label, 0.0im) for label in output_pairs]
    parameters = vcat(pack_complex(coefficients), [radius, spin, alpha])

    for potential in (:quartic, :sextic, :cosine)
        @testset "$potential" begin
            source = Dict(zip(output_pairs, unpack_complex(
                nonlinear_two_phase_source(parameters, source_basis; potential))))
            reference = nonlinear_source(single_cloud_values, radius, spin, alpha,
                polar_table, azimuthal_table, grid; modes=[(1, 1), (3, 3)], potential)

            @test isapprox(source[(1, (1, 0))], reference[(1, 1)]; atol=1e-13, rtol=1e-12)
            # A nonzero third harmonic is generated despite its absence from the input.
            @test abs(reference[(3, 3)]) > 1e-6
            @test isapprox(source[(3, (3, 0))], reference[(3, 3)]; atol=1e-13, rtol=1e-12)
            # This also has m=3, but its q_2=2 phase dependence cannot be generated.
            @test isapprox(source[(3, (-1, 2))], 0; atol=1e-13, rtol=0)
        end
    end
end
