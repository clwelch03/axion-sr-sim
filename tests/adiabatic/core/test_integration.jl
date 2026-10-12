using LinearAlgebra: norm

@testset "RK4 radial integration" begin
    # R'=Q, Q'=-R has the exact solution [cos(x), -sin(x)].
    function oscillator!(derivative, state, position)
        derivative[1] = state[2]
        derivative[2] = -state[1]
    end
    for (start, stop) in ((0.0, 1.0), (1.0, 0.0))
        initial = [cos(start), -sin(start)]
        saved_initial = copy(initial)
        exact = [cos(stop), -sin(stop)]
        errors = [norm(integrate_radial(oscillator!, initial, start, stop; step) - exact)
                  for step in (0.1, 0.05, 0.025)]
        # Fourth-order global error falls by 2^4 when the spacing is halved.
        @test 15 < errors[1] / errors[2] < 17
        @test 15 < errors[2] / errors[3] < 17
        @test initial == saved_initial
    end

    # A position-dependent RHS checks stage positions and a shortened spacing.
    # y'=3i*x^2 integrates exactly to y=i*x^3 under RK4.
    polynomial! = (derivative, state, position) -> (derivative .= 3im * position^2)
    for (start, stop) in ((0.0, 1.0), (1.0, 0.0))
        result = integrate_radial(polynomial!, [im * start^3], start, stop; step=0.13)
        @test isapprox(only(result), im * stop^3; atol=1e-14)
    end
    @test integrate_radial(oscillator!, [1, 0], 0.0, 0.0) == [1, 0]
    @test_throws ArgumentError integrate_radial(oscillator!, [1, 0], 0.0, 1.0; step=0)

    # Check every saved state against an analytic, position-dependent solution
    # in both directions. This also catches aliased snapshots of the final state.
    for (start, stop) in ((0.0, 1.0), (1.0, 0.0))
        profile = integrate_radial(polynomial!, [im * start^3], start, stop;
            step=0.13, save_profile=true)
        @test first(profile.positions) == start
        @test last(profile.positions) == stop
        @test vec(profile.states) ≈ im .* profile.positions .^ 3
        @test profile.states[:, end] == integrate_radial(
            polynomial!, [im * start^3], start, stop; step=0.13)
    end
end
