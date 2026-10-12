@testset "Coupled first-order state and harmonic frequencies" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    # Even m and an intentionally mixed ordering check m=n*m0, not m=n.
    modes = [(4, 4), (2, 2), (4, 2)]
    state = ComplexF64[0.1im, 0.2 + 0.1im, -0.3, 0.04, 0.05im, -0.06]
    saved_state = copy(state)
    derivative = similar(state)
    linear_tortoise_rhs!(derivative, state, 5.0, 0.5, 0.42, 0.3 + 0.01im,
                         modes, couplings; fundamental_m=2)
    @test state == saved_state
    @test derivative[1:3] == state[4:6]
    geometry = tortoise_geometry(5.0, 0.5)
    for (index, (l, m)) in enumerate(modes)
        radial = radial_operator_coefficients(geometry.radius, 0.5, 0.42, (m / 2) * (0.3 + 0.01im), l, m)
        mixed = sum(couplings[(l, other_l, m)] * state[j]
                    for (j, (other_l, other_m)) in enumerate(modes) if m == other_m)
        h = geometry.dr_drstar
        h_prime = (radial.delta_prime - 2geometry.radius * h) / geometry.radial_factor
        expected = h^2 * radial_second_derivative(state[index], state[3 + index] / h,
                                                  mixed, 0, radial) + h_prime * state[3 + index]
        @test isapprox(derivative[3 + index], expected; rtol=1e-12)
    end
end
