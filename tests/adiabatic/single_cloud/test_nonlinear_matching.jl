@testset "Nonlinear solve: fixed amplitude and quadratic frequency shift" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(5, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(5, polar_values, grid)
    spin, alpha = 0.99, 0.42
    modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
    linear = solve_linear_mode(0.4088 + 1e-7im, spin, alpha, [1, 3, 5], 1, couplings)
    scaled_horizon_seed = vcat(linear.horizon_amplitudes, zeros(ComplexF64, 3))
    scaled_outer_seed = vcat(linear.outer_amplitudes, zeros(ComplexF64, 3))
    real_shifts = Float64[]

    # Solve independently from the same linear branch at both amplitudes.
    # Dividing the matching residual by A0 prevents false convergence to zero.
    for amplitude in (1e-4, 5e-5)
        solution = solve_nonlinear_mode(
            linear.frequency, amplitude, spin, alpha, modes, couplings,
            polar_values, azimuthal_values, grid;
            horizon_guess=amplitude .* scaled_horizon_seed,
            outer_tail_guess=amplitude .* scaled_outer_seed[2:end],
            tolerance=1e-11,
        )
        @test solution.outer_amplitudes[1] == amplitude
        @test solution.scaled_residual_norm < 1e-11
        @test solution.iterations > 0
        @test isapprox(imag(solution.frequency), imag(linear.frequency); rtol=1e-3)

        # Re-evaluate the physical equations using the returned amplitudes.
        checked_residual = nonlinear_matching_residual(
            solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
            spin, alpha, modes, couplings, polar_values, azimuthal_values, grid,
        )
        @test norm(checked_residual) / amplitude < 1e-11
        @test solution.residual_norm ≈ amplitude * solution.scaled_residual_norm
        push!(real_shifts, real(solution.frequency - linear.frequency))
        println("A0=$amplitude: frequency=$(solution.frequency), " *
                "scaled residual=$(solution.scaled_residual_norm), " *
                "iterations=$(solution.iterations)")
        flush(stdout)
    end

    # The leading sine remainder is cubic in the field; its frequency correction
    # is quadratic. This checks weak-field scaling, not endpoint convergence.
    @test all(shift -> -1e-7 < shift < 0, real_shifts)
    @test isapprox(real_shifts[1], 4real_shifts[2]; rtol=2e-3)
    println("Real frequency shifts: $real_shifts; ratio=$(real_shifts[1] / real_shifts[2])")
end
