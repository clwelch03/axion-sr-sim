@testset "Fixed-scale two-phase profile coordinates" begin
    cloud_m = [1, 2]
    components = [fourier_component(q, [l], cloud_m) for (l,q) in
        [(0,(2,-1)), (2,(0,1)), (3,(-1,2)), (1,(1,0))]]
    amplitudes = [1e-4, 1e-5]
    scale = minimum(amplitudes)
    state = ComplexF64[0.408+1e-7im, 0.415+1e-9im, 1+2im, 3+4im, 5+6im, 7+8im, 9+10im, 11+12im]
    packed_state = pack_complex(state)
    grid = angular_grid(4, 5)
    polar_values = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar_values, grid)
    (; free_outer_coefficient_indices) = two_phase_coefficient_indices(components)
    settings = (; spin=0.99, alpha=0.42, components, cloud_m, couplings,
        polar_values, grid, phase1_nodes=grid.azimuth_nodes,
        phase2_nodes=grid.azimuth_nodes, potential=:quartic,
        coefficient_scale=scale, free_outer_coefficient_indices,
        initial_bound=BitVector([real(sum(component.q .* state[1:2])) < 0.42
            for component in components]),
        integration=(; inner=-30.0, outer=80.0, match=5.0, step=0.4,
            horizon_order=3, outer_order=7))
    problem = TwoPhaseProblem(; settings...)
    original = (copy(packed_state), copy(amplitudes))
    profile = unpack_two_phase_profile(packed_state, amplitudes, problem)
    @test profile isa TwoPhaseProfile
    @test profile.frequencies == state[1:2]
    @test profile.horizon_coefficients == scale .* state[3:6]
    @test profile.outer_coefficients[[4,2]] == amplitudes
    @test profile.outer_coefficients[[1,3]] == scale .* state[7:8]
    @test (packed_state, amplitudes) == original
    @test pack_two_phase_profile(profile, problem) ≈ packed_state

    # Repack in a different fixed scale; physical coefficients must not change.
    rescaled_problem = TwoPhaseProblem(; settings..., coefficient_scale=7scale)
    rescaled_state = pack_complex(vcat(state[1:2], state[3:end] ./ 7))
    rescaled_profile = unpack_two_phase_profile(rescaled_state, amplitudes, rescaled_problem)
    @test rescaled_profile.frequencies == profile.frequencies
    @test rescaled_profile.horizon_coefficients ≈ profile.horizon_coefficients
    @test rescaled_profile.outer_coefficients ≈ profile.outer_coefficients

    replacement_amplitudes = 2 .* amplitudes
    replacement_profile = unpack_two_phase_profile(packed_state, replacement_amplitudes, problem)
    @test replacement_profile.outer_coefficients[[4,2]] == replacement_amplitudes
    @test replacement_profile.outer_coefficients[[1,3]] == profile.outer_coefficients[[1,3]]
    @test replacement_profile.horizon_coefficients == profile.horizon_coefficients
    @test replacement_profile.frequencies == profile.frequencies

    function physical_coordinates(packed, norms)
        trial_profile = unpack_two_phase_profile(packed, norms, problem)
        return pack_complex(vcat(trial_profile.frequencies,
            trial_profile.horizon_coefficients, trial_profile.outer_coefficients))
    end
    amplitude_jacobian = ForwardDiff.jacobian(
        norms -> physical_coordinates(packed_state, norms), amplitudes)
    expected_amplitude_jacobian = zeros(20, 2)
    expected_amplitude_jacobian[10,1] = 1
    expected_amplitude_jacobian[8,2] = 1
    @test amplitude_jacobian == expected_amplitude_jacobian

    # Eight complex unknowns map into ten complex physical entries. The two
    # prescribed outer coefficients have zero derivatives with respect to state.
    expected_state_jacobian = zeros(20, 16)
    physical_rows = [1, 2, 3, 4, 5, 6, 7, 9]
    for (column, row) in enumerate(physical_rows)
        derivative = column <= 2 ? 1.0 : scale
        expected_state_jacobian[row, column] = derivative
        expected_state_jacobian[row + 10, column + 8] = derivative
    end
    @test ForwardDiff.jacobian(values -> physical_coordinates(values, amplitudes),
        packed_state) == expected_state_jacobian

    @test_throws DimensionMismatch unpack_two_phase_profile(packed_state[1:14], amplitudes, problem)
    @test_throws DimensionMismatch unpack_two_phase_profile(packed_state[1:15], amplitudes, problem)
    @test_throws DimensionMismatch unpack_two_phase_profile(packed_state, amplitudes[1:1], problem)
    @test_throws ArgumentError unpack_two_phase_profile(packed_state, [0.0, 1.0], problem)
    @test_throws ArgumentError unpack_two_phase_profile(packed_state, [Inf, 1.0], problem)
    @test_throws ArgumentError unpack_two_phase_profile(packed_state, amplitudes,
        TwoPhaseProblem(; settings..., coefficient_scale=0.0))
    @test_throws ArgumentError unpack_two_phase_profile(packed_state, amplitudes,
        TwoPhaseProblem(; settings..., coefficient_scale=NaN))
    @test_throws ArgumentError unpack_two_phase_profile(packed_state, amplitudes,
        TwoPhaseProblem(; settings..., components=[components[1:3]; components[1]]))
end
