@testset "Amplitude RK4 exact solutions and convergence" begin
    # Constant rate integrates exactly.
    @test rk4_amplitude_step(amplitude -> -0.5, 1.0, 0.2) ≈ 0.9

    # An amplitude initially at equilibrium stays there.
    @test rk4_amplitude_step(amplitude -> 2-amplitude, 2.0, 0.2) == 2.0

    function evolve_test(rate, step_count)
        amplitude = 1.0
        for _ in 1:step_count
            amplitude = rk4_amplitude_step(
                rate, amplitude, 1.0/step_count)
        end
        return amplitude
    end

    # Evolve to t=1 from A(0)=1:
    # dA/dt=A     gives A(t)=exp(t).
    # dA/dt=2-A   gives A(t)=2-exp(-t).
    for (rate, exact) in ((amplitude -> amplitude, exp(1.0)),
                          (amplitude -> 2-amplitude, 2-exp(-1.0)))
        coarse_error = abs(evolve_test(rate, 8) - exact)
        fine_error = abs(evolve_test(rate, 16) - exact)

        # Fourth-order global error: halving dt reduces error by ~16.
        @test 14 < coarse_error/fine_error < 18
        @test fine_error < 1e-6
    end

    @test_throws ArgumentError rk4_amplitude_step(identity, 1.0, 0.0)
    @test_throws ArgumentError rk4_amplitude_step(identity, 1.0, -0.1)
    @test_throws ArgumentError rk4_amplitude_step(identity, 1.0, Inf)
    @test_throws ErrorException rk4_amplitude_step(
        amplitude -> NaN, 1.0, 0.1)
end