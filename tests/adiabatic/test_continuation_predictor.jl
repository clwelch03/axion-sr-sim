@testset "Continuation prediction on unequal amplitude steps" begin
    # An affine complex branch is extrapolated exactly, including phase drift.
    exact(amplitude) = (
        frequency=0.4 + 1e-7im + amplitude * (-0.01 + 2e-8im),
        horizon_amplitudes=ComplexF64[1+2im, -3im] .+
            amplitude .* ComplexF64[0.2-0.3im, 0.7+0.1im],
        outer_amplitudes=ComplexF64[amplitude, 0.01im + amplitude*(0.2-0.1im)],
    )
    continuation = [(amplitude=amplitude, solution=exact(amplitude))
                    for amplitude in (0.1, 0.16)]
    unchanged = deepcopy(continuation)
    prediction = predict_nonlinear_mode(continuation, 0.19)
    @test prediction.frequency ≈ exact(0.19).frequency
    @test prediction.horizon_amplitudes ≈ exact(0.19).horizon_amplitudes
    @test prediction.outer_amplitudes ≈ exact(0.19).outer_amplitudes
    @test prediction.outer_amplitudes[1] == 0.19
    @test continuation == unchanged

    prediction = predict_nonlinear_mode(continuation[1:1], 0.2)
    @test prediction.frequency == continuation[1].solution.frequency
    @test prediction.horizon_amplitudes ≈ 2continuation[1].solution.horizon_amplitudes
    @test prediction.outer_amplitudes ≈ 2continuation[1].solution.outer_amplitudes
    @test_throws ArgumentError predict_nonlinear_mode(continuation, 0.16)
    @test_throws ArgumentError predict_nonlinear_mode(reverse(continuation), 0.2)
    @test_throws ArgumentError predict_nonlinear_mode(continuation, Inf)
    @test_throws ArgumentError predict_nonlinear_mode([], 0.2)
end

@testset "Cubic continuation prediction and short-history fallback" begin
    exact(amplitude) = (
        frequency=0.4 - 0.01amplitude^2 + 1e-7im * (1 + amplitude^3),
        horizon_amplitudes=ComplexF64[1 + im*amplitude^3, amplitude^2 - 2im*amplitude],
        outer_amplitudes=ComplexF64[amplitude, (1+im)*amplitude^3 - amplitude^2],
    )
    continuation = [(; amplitude, solution=exact(amplitude)) for amplitude in (0.1, 0.2, 0.31, 0.45)]
    unchanged = deepcopy(continuation)
    prediction = predict_nonlinear_mode_cubic(continuation, 0.5; history_spacing=0.1)
    @test prediction.frequency ≈ exact(0.5).frequency
    @test prediction.horizon_amplitudes ≈ exact(0.5).horizon_amplitudes
    @test prediction.outer_amplitudes ≈ exact(0.5).outer_amplitudes
    @test prediction.outer_amplitudes[1] == 0.5
    @test continuation == unchanged
    @test predict_nonlinear_mode_cubic(continuation[1:2], 0.3) ==
        predict_nonlinear_mode(continuation[1:2], 0.3)
    # Fixed target offsets select 0.26 twice despite sufficient earlier history.
    uneven = [(; amplitude, solution=exact(amplitude)) for amplitude in (0.1, 0.2, 0.26, 0.35, 0.36)]
    prediction = predict_nonlinear_mode_cubic(uneven, 0.4; history_spacing=0.05)
    @test prediction.frequency ≈ exact(0.4).frequency
    @test prediction.horizon_amplitudes ≈ exact(0.4).horizon_amplitudes
    @test prediction.outer_amplitudes ≈ exact(0.4).outer_amplitudes
    @test prediction != predict_nonlinear_mode(uneven, 0.4)
    @test_throws ArgumentError predict_nonlinear_mode_cubic(continuation, 0.5; history_spacing=0)
end
