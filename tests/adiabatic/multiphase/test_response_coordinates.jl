@testset "Fixed-scale two-phase trial unpacking" begin
    cloud_m = [1, 2]
    components = [fourier_component(q, [l], cloud_m) for (l,q) in
        [(0,(2,-1)), (2,(0,1)), (3,(-1,2)), (1,(1,0))]]
    amplitudes = [1e-4, 1e-5]
    scale = minimum(amplitudes)
    state = ComplexF64[0.408+1e-7im, 0.415+1e-9im, 1+2im, 3+4im, 5+6im, 7+8im, 9+10im, 11+12im]
    original = (copy(state), copy(amplitudes))
    frequencies, horizon, outer = unpack_two_phase_trial(state, amplitudes, components; coefficient_scale=scale)
    @test frequencies == state[1:2]
    @test horizon == scale .* state[3:6]
    @test outer[[4,2]] == amplitudes
    @test outer[[1,3]] == scale .* state[7:8]
    @test (state, amplitudes) == original
    # Repack in a different fixed scale; physical coefficients must not change.
    rescaled_state = vcat(state[1:2], state[3:end] ./ 7)
    @test all(isapprox.(unpack_two_phase_trial(rescaled_state, amplitudes, components;
        coefficient_scale=7scale), (frequencies, horizon, outer)))
    physical = (packed, norms) -> pack_complex(vcat(unpack_two_phase_trial(
        unpack_complex(packed), norms, components; coefficient_scale=scale)...))
    packed = pack_complex(state)
    amplitude_jacobian = ForwardDiff.jacobian(norms -> physical(packed, norms), amplitudes)
    expected = zeros(20, 2)
    expected[10,1] = 1
    expected[8,2] = 1
    @test amplitude_jacobian == expected
    @test all(isfinite, ForwardDiff.jacobian(values -> physical(values, amplitudes), packed))
    @test_throws DimensionMismatch unpack_two_phase_trial(state[1:7], amplitudes, components; coefficient_scale=scale)
    @test_throws DimensionMismatch unpack_two_phase_trial(state, amplitudes[1:1], components; coefficient_scale=scale)
    @test_throws ArgumentError unpack_two_phase_trial(state, [0.0, 1.0], components; coefficient_scale=scale)
    @test_throws ArgumentError unpack_two_phase_trial(state, amplitudes, components; coefficient_scale=0.0)
    @test_throws ArgumentError unpack_two_phase_trial(state, amplitudes, [components[1:3]; components[1]]; coefficient_scale=scale)
end
