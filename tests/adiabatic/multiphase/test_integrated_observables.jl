# Candidate checkpoint: load observables.jl and the radial wrapper before this file.
# Register after user insertion and loader inclusion.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "Radial observables: exact cubic energy and branch orientation" begin
    grid = angular_grid(8, 9)
    polar = build_polar_values(1, grid.cosine_theta_nodes)
    cloud_m, frequencies = [1, 2], [0.4, 0.41]
    components = [fourier_component((1, 0), [1], cloud_m)]
    spin = 0.7
    function branch(start, stop, intervals)
        positions = collect(range(start, stop; length=intervals+1))
        # Manufactured snapshots, not an ODE solution: R=0 and prescribed Q
        # give exactly dE/dr_*=1+r_*^3, J=0 and both currents=0.
        states = zeros(ComplexF64, 2, length(positions))
        for (index, position) in enumerate(positions)
            geometry = tortoise_geometry(position, spin)
            states[2, index] = sqrt((1+position^3)/(2pi*geometry.radial_factor))
        end
        return (; positions, states)
    end
    expected = (8.0+8.0^4/4) - (3.0+3.0^4/4)
    for (inner_count, outer_count) in ((2, 3), (3, 4), (5, 6)), reverse_branches in (false, true)
        inner = reverse_branches ? branch(5, 3, inner_count) : branch(3, 5, inner_count)
        outer = reverse_branches ? branch(5, 8, outer_count) : branch(8, 5, outer_count)
        result = two_phase_integrated_observables((; from_horizon=inner, from_outer=outer),
            spin, 0.42, frequencies, components, cloud_m, polar, grid,
            grid.azimuth_nodes, [0.0])
        @test isapprox(result.energy, expected; rtol=2e-12)
        # Raising the tensor index cancels two Kerr off-diagonal contributions;
        # test zero angular momentum to roundoff relative to the charge scale.
        @test abs(result.angular_momentum) < 10eps(Float64)*expected
        @test all(iszero, (result.horizon_energy_flux,
            result.outer_energy_flux, result.horizon_angular_momentum_flux,
            result.outer_angular_momentum_flux))
    end
    println("Integrated manufactured cubic energy: expected ", expected)
end

@testset "Radial observables: angular momentum, endpoint signs, and old energy" begin
    grid = angular_grid(12, 17)
    polar = build_polar_values(1, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(1, grid.azimuth_nodes)
    cloud_m, frequencies = [1, 2], [0.4, 0.41]
    components = [fourier_component((1, 0), [1], cloud_m)]
    spin, alpha, wave_number = 0.0, 0.42, 0.07
    function branch(start, stop, intervals)
        positions = collect(range(start, stop; length=intervals+1))
        states = zeros(ComplexF64, 2, length(positions))
        for (index, position) in enumerate(positions)
            radius = tortoise_geometry(position, spin).radius
            # Prescribed snapshots: exact J density 1+r_*^3, positive outward flux.
            states[1, index] = sqrt((1+position^3)/(4pi*frequencies[1]*radius^2))
            states[2, index] = im * wave_number * states[1, index]
        end
        return (; positions, states)
    end
    profiles = (; from_horizon=branch(3, 5, 4), from_outer=branch(8, 5, 6))
    expected_momentum = (8.0+8.0^4/4) - (3.0+3.0^4/4)
    for potential in (:cosine, :quartic, :sextic)
        result = two_phase_integrated_observables(profiles, spin, alpha, frequencies,
            components, cloud_m, polar, grid, -grid.azimuth_nodes, [0.0]; potential)
        expected_energy = sum(integrate_energy_branch(branch, frequencies[1], spin,
            alpha, [(1, 1)], polar, azimuthal, grid; potential)
            for branch in (profiles.from_horizon, profiles.from_outer))
        @test isapprox(result.energy, expected_energy; rtol=2e-12)
        @test isapprox(result.angular_momentum, expected_momentum; rtol=2e-12)
        @test isapprox(result.horizon_energy_flux, -wave_number*(1+3^3); rtol=2e-12)
        @test isapprox(result.outer_energy_flux, wave_number*(1+8^3); rtol=2e-12)
        @test isapprox(result.horizon_angular_momentum_flux,
            -wave_number*(1+3^3)/frequencies[1]; rtol=2e-12)
        @test isapprox(result.outer_angular_momentum_flux,
            wave_number*(1+8^3)/frequencies[1]; rtol=2e-12)
    end
end

@testset "Radial observables: saved matching histories and validation" begin
    grid = angular_grid(8, 9)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar, grid)
    cloud_m, frequencies = [1, 2], [0.408, 0.416]
    components = [fourier_component(q, [l], cloud_m) for (l, q) in
        [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]]
    spin, alpha = 0.7, 0.42
    # At these short endpoints the driven boundary series strongly amplifies
    # coefficients. Keep the trial dilute; this checks the saved-history API,
    # not boundary-series accuracy or physical convergence.
    profiles = two_phase_matching_residual(frequencies, fill(1e-13+2e-14im, 4),
        fill(1e-12-3e-13im, 4), spin, alpha, components, cloud_m, couplings,
        polar, grid, grid.azimuth_nodes, grid.azimuth_nodes;
        inner=3.0, match=4.0, outer=5.0, step=0.5, save_profile=true, potential=:quartic)
    # These are trial histories, not a converged physical profile.
    evaluate = (histories, cs=components) -> two_phase_integrated_observables(histories,
        spin, alpha, frequencies, cs, cloud_m, polar, grid,
        grid.azimuth_nodes, grid.azimuth_nodes; potential=:quartic)
    original = deepcopy(profiles)
    reference = evaluate(profiles)
    @test all(isfinite, values(reference))
    reverse_branch = branch -> (; positions=reverse(branch.positions), states=branch.states[:, end:-1:1])
    reversed = evaluate((; from_horizon=reverse_branch(profiles.from_horizon),
        from_outer=reverse_branch(profiles.from_outer)))
    @test collect(values(reversed)) ≈ collect(values(reference))
    permutation = [3, 1, 4, 2]
    rows = vcat(permutation, 4 .+ permutation)
    reorder = branch -> (; branch.positions, states=branch.states[rows, :])
    reordered = evaluate((; from_horizon=reorder(profiles.from_horizon),
        from_outer=reorder(profiles.from_outer)), components[permutation])
    @test collect(values(reordered)) ≈ collect(values(reference))
    @test profiles == original

    bad_positions = positions -> merge(profiles,
        (; from_horizon=(; positions, states=profiles.from_horizon.states)))
    @test_throws ArgumentError evaluate(bad_positions([3.0, 3.6, 4.0]))
    @test_throws ArgumentError evaluate(bad_positions([3.0, 4.0, 3.5]))
    @test_throws ArgumentError evaluate(bad_positions([3.0, 3.0, 4.0]))
    @test_throws ArgumentError evaluate(bad_positions([3.0, NaN, 4.0]))
    @test_throws ArgumentError evaluate(bad_positions([3.0, 3.4, 3.8]))
    @test_throws ArgumentError evaluate(merge(profiles, (; from_horizon=(;
        positions=[3.0, 4.0], states=profiles.from_horizon.states[:, 1:2]))))
    @test_throws DimensionMismatch evaluate(merge(profiles, (; from_horizon=(;
        profiles.from_horizon.positions, states=profiles.from_horizon.states[1:7, :]))))
    nonfinite_states = copy(profiles.from_horizon.states)
    nonfinite_states[1, 2] = NaN
    @test_throws ArgumentError evaluate(merge(profiles, (; from_horizon=(;
        profiles.from_horizon.positions, states=nonfinite_states))))
    overflowing_states = copy(profiles.from_horizon.states)
    overflowing_states[1, 2] = 1e200
    @test_throws ArgumentError evaluate(merge(profiles, (; from_horizon=(;
        profiles.from_horizon.positions, states=overflowing_states))))
end
