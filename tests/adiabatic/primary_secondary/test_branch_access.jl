using Serialization

@testset "Shared primary/secondary branch access" begin
    # Distinct row tolerances catch accidental use of the branch default.
    solution = (frequency=0.298 + 1e-12im,
        horizon_amplitudes=ComplexF64[0.2, 0.01],
        outer_amplitudes=ComplexF64[1, 0.03],
        residual=ComplexF64[1e-7, 0, 0, 0], residual_norm=1e-7, iterations=4)
    rows = [
        (; primary_index=12, amplitude=0.02, omega1=0.296 + 1e-8im,
           solution, primary_residual_change=2e-13, tolerance=2e-6),
        (; primary_index=19, amplitude=0.04, omega1=0.295 + 1e-8im,
           solution, primary_residual_change=3e-13, tolerance=5e-6),
    ]
    legacy_secondary = (; primary_path="original-primary.jls",
        spin=0.99, alpha=0.3, secondary_modes=[(2, 2), (4, 2)],
        primary_options=(; step=0.1, inner=-100.0, outer=300.0, match=5.0),
        tolerance=1e-8, reading="averaged", linear=solution, rows)
    nonlinear_solution = (; solution..., scaled_residual_norm=5e-6)
    legacy_primary = (; spin=0.99, alpha=0.3, potential="cosine",
        continuation=[(; amplitude=0.02, solution=nonlinear_solution)],
        rows=[(; amplitude=0.02, energy=600.0)], tolerance=1e-5)

    mktempdir() do directory
        primary_path, secondary_path = joinpath.(directory, ("primary.jls", "secondary.jls"))
        serialize(primary_path, legacy_primary)
        serialize(secondary_path, legacy_secondary)
        primary = read_branch(primary_path;
            checkpoint_adapter=single_cloud_continuation_checkpoint_adapter)
        secondary = read_branch(secondary_path;
            checkpoint_adapter=secondary_continuation_checkpoint_adapter)

        @test primary isa ContinuationBranch
        @test secondary isa ContinuationBranch
        @test branch_point(primary, 1).diagnostics.observables.energy == 600
        @test branch_point(secondary, 1).parameter_values.amplitude == 1
        @test branch_point(secondary, 1).parameter_values.primary_amplitude == 0.02
        @test branch_point(secondary, 1).backgrounds.primary.point_id == 12
        @test branch_point(secondary, 2).backgrounds.primary.point_id == 19
        @test branch_point(secondary, 2).backgrounds.primary.checkpoint_path == "original-primary.jls"
        @test branch_point(secondary, 2).diagnostics.tolerance == 5e-6
        @test secondary.settings.linear == solution

        # IDs continue to select the same physical calculation after reordering.
        reverse!(secondary.points)
        @test branch_point(secondary, 1).diagnostics.tolerance == 2e-6
        @test_throws ArgumentError branch_point(secondary, 100)
        duplicate = ContinuationBranch(secondary.settings,
            [secondary.points[1], secondary.points[1]], secondary.provenance)
        @test_throws ArgumentError branch_point(duplicate, 2)

        # A native branch can be read without an adapter; coefficients survive IO.
        native_path = joinpath(directory, "native.jls")
        serialize(native_path, secondary)
        restored = read_branch(native_path)
        @test branch_point(restored, 2).profile.outer_coefficients == solution.outer_amplitudes
        @test branch_point(restored, 2).diagnostics.residual == solution.residual
        @test branch_point(restored, 2).diagnostics.primary_frequency == rows[2].omega1
        @test_throws ArgumentError read_branch(secondary_path)
    end

    # Older files without per-point tolerances must retain that uncertainty.
    older_row = (; (name => value for (name, value) in pairs(rows[1]) if name != :tolerance)...)
    older = secondary_continuation_checkpoint_adapter(
        (; legacy_secondary..., rows=[older_row]); checkpoint_path="older.jls")
    @test only(older.points).diagnostics.tolerance === nothing
    only(older.points).profile.outer_coefficients[2] = 99
    @test solution.outer_amplitudes[2] == 0.03

    nonunit_solution = (; solution..., outer_amplitudes=ComplexF64[2, 0.03])
    nonunit = (; legacy_secondary..., rows=[(; rows[1]..., solution=nonunit_solution)])
    @test_throws ArgumentError secondary_continuation_checkpoint_adapter(nonunit; checkpoint_path="bad.jls")

    empty_branch = secondary_continuation_checkpoint_adapter(
        (; legacy_secondary..., rows=[]); checkpoint_path="empty.jls")
    @test isempty(empty_branch.points)
    @test empty_branch.settings.linear == solution
end

@testset "Linear matching with retained radial samples" begin
    grid = angular_grid(12, 25)
    polar_values = build_polar_values(4, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(4, polar_values, grid)
    # Exercise both fundamental-m conventions without an expensive eigenmode solve.
    for (m, frequency) in ((1, 0.296 + 1e-8im), (2, 0.298 + 1e-12im))
        l_values = [m, m + 2]
        horizon_coefficients, outer_coefficients = ComplexF64[0.1, 0.01], ComplexF64[1, 0.03]
        options = (; step=1.0, inner=-20.0, outer=50.0, match=5.0)
        mismatch = matching_residual(frequency, horizon_coefficients, outer_coefficients,
            0.99, 0.3, l_values, m, couplings; options...)
        radial_result = matching_residual(frequency, horizon_coefficients, outer_coefficients,
            0.99, 0.3, l_values, m, couplings; options..., save_profile=true)
        @test radial_result.residual == mismatch
        @test radial_result.from_horizon.positions[[1, end]] == [-20, 5]
        @test radial_result.from_outer.positions[[1, end]] == [50, 5]
        @test size(radial_result.from_horizon.states, 1) == 4
        @test all(isfinite, radial_result.from_outer.states)
    end
end
