# Candidate checkpoint: load two_phase_potential_energy_density before this file.
# Register after the user inserts the tested helper into observables.jl.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "Joint potential energy at one radius" begin
    grid = angular_grid(24, 32)
    polar = build_polar_values(2, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(2, grid.azimuth_nodes)
    phase1_nodes = -grid.azimuth_nodes
    phase2_nodes = 2pi .* (0:34) ./ 35
    position, spin, alpha = 8.0, 0.7, 0.42
    first_value, second_value = 0.8+0.3im, 0.4-0.2im
    first_values = Dict((1, (1, 0)) => first_value)
    second_values = Dict((2, (0, 1)) => second_value)
    joint_values = merge(first_values, second_values)
    density = (values, potential) -> two_phase_potential_energy_density(
        values, position, spin, alpha, [1, 2], polar, grid,
        phase1_nodes, phase2_nodes; potential)

    for potential in (:cosine, :quartic, :sextic)
        state = ComplexF64[first_value, 0]
        reference = energy_per_tortoise_radius(state, position, 0.4, spin,
            alpha, [(1, 1)], polar, azimuthal, grid; potential) -
            energy_per_tortoise_radius(state, position, 0.4, spin,
                0.0, [(1, 1)], polar, azimuthal, grid; potential)
        @test isapprox(density(first_values, potential), reference; rtol=1e-12)

        shifted_values = Dict((l, q) => value * cis(-q[1]*0.37 + q[2]*0.53)
            for ((l, q), value) in joint_values)
        @test isapprox(density(shifted_values, potential),
            density(joint_values, potential); rtol=1e-12)
    end

    # For U=Phi^2/2-Phi^4/24, the independent-phase interaction is
    # <U(Phi1+Phi2)-U(Phi1)-U(Phi2)> = -|f1|^2 |f2|^2.
    # With f1=u*Y11 and f2=v*Y22, the exact polar integral is
    # |u*v|^2 * (9r^2+spin^2)/14, including rho^2.
    geometry = tortoise_geometry(position, spin)
    expected_interaction = -2pi * geometry.dr_drstar * alpha^2 *
        abs2(first_value) * abs2(second_value) *
        (9geometry.radius^2 + spin^2)/14
    interaction = density(joint_values, :quartic) -
        density(first_values, :quartic) - density(second_values, :quartic)
    @test isapprox(interaction, expected_interaction; rtol=1e-12)
    @test iszero(density(Dict(), :cosine))
    @test iszero(two_phase_potential_energy_density(joint_values, position,
        spin, 0.0, [1, 2], polar, grid, phase1_nodes, phase2_nodes))
    @test_throws ArgumentError two_phase_potential_energy_density(joint_values,
        position, spin, alpha, [1, 2], polar, grid, Float64[], phase2_nodes)
    println("Quartic interaction density: ", interaction,
        "; exact: ", expected_interaction)
end
