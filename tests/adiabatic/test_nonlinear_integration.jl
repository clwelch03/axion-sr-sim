@testset "Nonlinear matching: linear limit and radial convergence" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(3, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(3, polar_values, grid)
    spin, alpha, frequency = 0.99, 0.42, 0.4088 + 1.5e-7im
    modes = [(1, 1), (3, 3)]
    horizon_amplitudes = ComplexF64[0.02, 0]
    outer_amplitudes = ComplexF64[0.015im, 0]
    saved_horizon, saved_outer = copy(horizon_amplitudes), copy(outer_amplitudes)

    # Short intervals and leading boundary terms isolate the integration check;
    # these endpoints are not intended to approximate a physical eigenmode.
    function mismatch(scale; step=0.1)
        return nonlinear_matching_residual(
            frequency, scale * horizon_amplitudes, scale * outer_amplitudes,
            spin, alpha, modes, couplings, polar_values, azimuthal_values, grid;
            inner=-10.0, outer=20.0, match=5.0, horizon_order=0, outer_order=0, step,
        )
    end

    linear_rhs! = (derivative, state, position) -> linear_tortoise_rhs!(
        derivative, state, position, spin, alpha, frequency, modes, couplings)
    inner_state = boundary_state(:horizon, -10.0, spin, alpha, frequency,
        modes, couplings, horizon_amplitudes; horizon_order=0)
    outer_state = boundary_state(:outer, 20.0, spin, alpha, frequency,
        modes, couplings, outer_amplitudes; outer_order=0)
    linear_mismatch = integrate_radial(linear_rhs!, inner_state, -10.0, 5.0; step=0.1) -
        integrate_radial(linear_rhs!, outer_state, 20.0, 5.0; step=0.1)
    @test isapprox(mismatch(1e-5) / 1e-5, linear_mismatch; rtol=1e-8)
    @test mismatch(0.0) == zeros(ComplexF64, 4)

    # A cubic source creates m=3 despite zero m=3 endpoint amplitudes.
    # Halving all input amplitudes must reduce this signal by approximately eight.
    weak = mismatch(0.01)
    half_weak = mismatch(0.005)
    @test linear_mismatch[[2, 4]] == zeros(ComplexF64, 2)
    @test norm(weak[[2, 4]]) > 1e-12
    @test isapprox(weak[[2, 4]], 8half_weak[[2, 4]]; rtol=1e-3)

    coarse, medium, fine = (mismatch(1.0; step) for step in (0.4, 0.2, 0.1))
    convergence_ratio = norm(coarse - medium) / norm(medium - fine)
    @test 12 < convergence_ratio < 20
    @test all(isfinite, fine)
    @test horizon_amplitudes == saved_horizon
    @test outer_amplitudes == saved_outer
    profile = nonlinear_matching_residual(
        frequency, horizon_amplitudes, outer_amplitudes,
        spin, alpha, modes, couplings, polar_values, azimuthal_values, grid;
        inner=-10.0, outer=20.0, match=5.0, horizon_order=0, outer_order=0,
        step=0.1, save_profile=true,
    )
    @test profile.residual == fine
    @test profile.from_horizon.states[:, 1] == inner_state
    @test profile.from_outer.states[:, 1] == outer_state
    @test last(profile.from_horizon.positions) == last(profile.from_outer.positions) == 5.0
    @test_throws ArgumentError nonlinear_matching_residual(
        frequency, horizon_amplitudes, outer_amplitudes,
        spin, alpha, modes, couplings, polar_values, azimuthal_values, grid;
        match=100.0,
    )
    println("Nonlinear matching RK4 refinement ratio: $convergence_ratio")
end
