@testset "Two-phase amplitude prediction and normalization" begin
    labels = [(1, (1, 0)), (2, (0, 1)), (0, (2, -1)), (3, (-1, 2))]
    profile = TwoPhaseProfile(ComplexF64[0.4, 0.41],
        ComplexF64[1+im, 2-im, 3+2im, 4-im], ComplexF64[0.1, 0.01, 3im, 4im])
    seed = (; profile, cloud_amplitudes=[0.1, 0.01],
        residual_weights=collect(1.0:16.0),
        problem=(; source_basis=(; mode_labels=labels)))
    original = deepcopy(seed)
    prediction = predict_two_phase_amplitude(seed, 0.02)
    @test prediction.cloud_amplitudes == [0.1, 0.02]
    @test prediction.profile.frequencies == profile.frequencies
    @test prediction.profile.horizon_coefficients == [1, 2, 2, 4] .* profile.horizon_coefficients
    @test prediction.profile.outer_coefficients == [1, 2, 2, 4] .* profile.outer_coefficients
    @test prediction.residual_weights == repeat([1, 2, 2, 2], 4) .* seed.residual_weights
    @test seed.cloud_amplitudes == original.cloud_amplitudes
    @test seed.residual_weights == original.residual_weights
    @test profile.horizon_coefficients == original.profile.horizon_coefficients
    @test profile.outer_coefficients == original.profile.outer_coefficients
    @test prediction.profile.frequencies !== profile.frequencies
    primary = predict_two_phase_amplitude(seed, 0.2; cloud_index=1)
    @test primary.cloud_amplitudes == [0.2, 0.01]
    @test primary.profile.horizon_coefficients == [2, 1, 4, 2] .* profile.horizon_coefficients
    @test primary.residual_weights == repeat([2, 1, 2, 2], 4) .* seed.residual_weights
    @test_throws ArgumentError predict_two_phase_amplitude(seed, 0.02; cloud_index=3)
    @test_throws ArgumentError predict_two_phase_amplitude(seed, Inf)
end
