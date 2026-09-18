@testset "Nonlinear radial RHS: sign, state ordering, and linear limit" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    spin, alpha, frequency, radius = 0.99, 0.42, 0.4088 + 1.5e-7im, 4.0
    coordinate = tortoise_coordinate(radius - kerr_horizons(spin).r_plus, spin)
    # Deliberately put the generated harmonic first, rather than the fundamental.
    modes = [(3, 3), (1, 1), (3, 1)]
    amplitude = 1e-3 * (1 + 0.4im)
    state = ComplexF64[0, amplitude, 0, 0, 0.02amplitude, -0.01amplitude]
    saved_state = copy(state)

    function derivatives(state, position; mass=alpha, active_modes=modes, fundamental_m=1)
        linear_derivative, nonlinear_derivative = similar(state), similar(state)
        linear_tortoise_rhs!(linear_derivative, state, position, spin, mass, frequency,
                             active_modes, couplings; fundamental_m)
        nonlinear_tortoise_rhs!(nonlinear_derivative, state, position, spin, mass, frequency,
            active_modes, couplings, polar, azimuthal, grid; fundamental_m)
        return linear_derivative, nonlinear_derivative
    end

    linear, nonlinear = derivatives(state, coordinate)
    correction = nonlinear - linear
    @test state == saved_state
    @test nonlinear[1:3] == state[4:6]
    @test nonlinear[1:3] == linear[1:3]

    # Independent cubic predictions, including the minus sign and dr_* conversion.
    source_scale = (radius^2 - 2radius + spin^2) / (radius^2 + spin^2)^2
    source_11 = alpha^2 * abs2(amplitude) * amplitude * (3radius^2/10 + 3spin^2/70)
    source_33 = alpha^2 * amplitude^3 * sqrt(3/280) * (radius^2 + spin^2/9)
    @test isapprox(correction[5], -source_scale * source_11; rtol=2e-6, atol=0)
    @test isapprox(correction[4], -source_scale * source_33; rtol=2e-6, atol=0)
    @test linear[4] == 0  # No m=3 field or slope, hence no linear m=3 forcing.
    @test real(nonlinear[4]) < 0  # The nonlinear source generates that harmonic.

    half_linear, half_nonlinear = derivatives(state / 2, coordinate)
    half_correction = half_nonlinear - half_linear
    for index in (4, 5)
        @test isapprox(correction[index], 8half_correction[index]; rtol=2e-6, atol=0)
    end
    # Cubic forcing divided by the linear RHS falls quadratically in amplitude.
    relative_change = norm(correction) / norm(linear)
    half_relative_change = norm(half_correction) / norm(half_linear)
    @test isapprox(relative_change, 4half_relative_change; rtol=2e-6, atol=0)

    # Nonzero Q with zero R must not create a potential source.
    slopes_only = copy(state)
    slopes_only[1:3] .= 0
    zero_linear, zero_nonlinear = derivatives(slopes_only, coordinate)
    @test zero_nonlinear == zero_linear

    # mu=0 removes the potential entirely; also check fundamental_m forwarding.
    massless_linear, massless_nonlinear = derivatives(state, coordinate;
        mass=0.0, active_modes=[(4, 4), (2, 2), (4, 2)], fundamental_m=2)
    @test massless_nonlinear == massless_linear

    # Stable Delta must still suppress the source deep in the near-horizon region.
    _, near_horizon = derivatives(state, -400.0)
    geometry = tortoise_geometry(-400.0, spin)
    near_source_33 = alpha^2 * amplitude^3 * sqrt(3/280) *
                     (geometry.radius^2 + spin^2/9)
    expected = -geometry.delta / geometry.radial_factor^2 * near_source_33
    @test all(isfinite, near_horizon)
    @test isapprox(near_horizon[4], expected; rtol=2e-6, atol=0)
end
