# Check the library matching residual against the single-cloud implementation,
# and verify joint ordering and short radial integration in both directions.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "Two-phase matching single-cloud limit" begin
    grid = angular_grid(12, 24)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(3, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(3, polar, grid)
    cloud_m, frequencies = [1, 2], [0.405+1e-7im, 0.415+2e-7im]
    components = [fourier_component((1, 0), [1, 3], cloud_m),
                  fourier_component((3, 0), [3], cloud_m)]
    modes = [(1, 1), (3, 1), (3, 3)]
    horizon = ComplexF64[1e-4+2e-5im, -2e-6+1e-6im, 1e-6-1e-7im]
    outer = ComplexF64[1e-5, 2e-7+1e-7im, -1e-7+2e-7im]
    options = (; inner=-10.0, match=5.0, outer=30.0, step=0.4)
    for potential in (:quartic, :sextic, :cosine)
        result = two_phase_matching_residual(frequencies, horizon, outer,
            0.7, 0.42, components, cloud_m, couplings, polar, grid,
            -grid.azimuth_nodes, 2pi .* (0:4) ./ 5; potential, options...,
            save_profile=true)
        reference = nonlinear_matching_residual(frequencies[1], horizon, outer,
            0.7, 0.42, modes, couplings, polar, azimuthal, grid;
            potential, options..., save_profile=true)
        # Check each entire integration history, not just a cancelling difference.
        @test isapprox(result.from_horizon.states, reference.from_horizon.states;
            atol=1e-12, rtol=1e-11)
        @test isapprox(result.from_outer.states, reference.from_outer.states;
            atol=1e-12, rtol=1e-11)
        @test isapprox(result.residual, reference.residual; atol=1e-12, rtol=1e-11)
        plain = two_phase_matching_residual(frequencies, horizon, outer,
            0.7, 0.42, components, cloud_m, couplings, polar, grid,
            -grid.azimuth_nodes, 2pi .* (0:4) ./ 5; potential, options...)
        @test plain == result.residual
    end
end

@testset "Joint matching ordering and radial refinement" begin
    grid = angular_grid(12, 9)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar, grid)
    phase1, phase2 = grid.azimuth_nodes, 2pi .* (0:10) ./ 11
    cloud_m, frequencies = [1, 2], [0.405+1e-7im, 0.415+2e-7im]
    components = [fourier_component((1, 0), [1, 3], cloud_m),
                  fourier_component((0, 1), [2], cloud_m),
                  fourier_component((2, -1), [0], cloud_m),
                  fourier_component((-1, 2), [3], cloud_m)]
    horizon = ComplexF64[1e-4+2e-5im, -2e-6+1e-6im, 8e-5-1e-5im, 0, 0]
    outer = ComplexF64[1e-5, 2e-7+1e-7im, 2e-6, 0, 0]
    horizon_original, outer_original = copy(horizon), copy(outer)
    options = (; inner=-10.0, match=5.0, outer=30.0, step=0.4, potential=:quartic)
    evaluate = (cs, hs, bs; kwargs...) -> two_phase_matching_residual(
        frequencies, hs, bs, 0.7, 0.42, cs, cloud_m, couplings,
        polar, grid, phase1, phase2; kwargs...)
    result = evaluate(components, horizon, outer; options..., save_profile=true)
    @test all(isfinite, result.residual)
    @test length(result.residual) == 10
    @test horizon == horizon_original && outer == outer_original
    @test (first(result.from_horizon.positions), last(result.from_horizon.positions)) == (-10.0, 5.0)
    @test (first(result.from_outer.positions), last(result.from_outer.positions)) == (30.0, 5.0)
    # Unequal block sizes make an incorrect concatenation of [R; Q] blocks visible.
    reordered = [components[4], fourier_component((1, 0), [3, 1], cloud_m),
                 components[3], components[2]]
    permutation = [5, 2, 1, 4, 3]
    reordered_residual = evaluate(reordered, horizon[permutation], outer[permutation]; options...)
    @test isapprox(reordered_residual, result.residual[vcat(permutation, 5 .+ permutation)];
        atol=1e-12, rtol=1e-11)
    @test_throws ArgumentError evaluate(components, horizon, outer; inner=5.0, match=5.0, outer=30.0)
    @test_throws ArgumentError evaluate(components[1:0], ComplexF64[], ComplexF64[]; options...)
    @test_throws DimensionMismatch evaluate(components, horizon[1:4], outer; options...)
    @test_throws DimensionMismatch evaluate(components, horizon, outer[1:4]; options...)

    # Isolate RK4 convergence from asymptotic-boundary accuracy: integrate a
    # fixed, finite nonlinear state over a short interval, in both directions.
    initial = ComplexF64[0.2+0.1im, 0.02-0.01im, 0.15-0.07im, 0, 0,
                        0.01-0.02im, 0.003im, -0.02+0.01im, 0, 0]
    rhs! = (derivative, state, position) -> nonlinear_two_phase_rhs!(derivative,
        state, position, 0.7, 0.42, components, frequencies, cloud_m, couplings,
        polar, grid, phase1, phase2; potential=:quartic)
    for (start, stop) in ((4.0, 8.0), (8.0, 4.0))
        endpoints = [integrate_radial(rhs!, initial, start, stop; step)
            for step in (0.4, 0.2, 0.1)]
        differences = [maximum(abs.(endpoints[j] - endpoints[j+1])) for j in 1:2]
        ratio = differences[1] / differences[2]
        @test 12 < ratio < 20 # Fourth-order global error predicts a ratio of 16.
        println("Two-phase RK4 ", start, " -> ", stop,
            ": endpoint differences=", differences, ", ratio=", ratio)
    end
end
