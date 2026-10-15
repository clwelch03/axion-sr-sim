# Candidate checkpoint: load kerr_metric, observables.jl, and
# two_phase_observables_at_radius before this file. Register after insertion.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))
using ForwardDiff

@testset "Combined observables: single-cloud energy" begin
    grid = angular_grid(24, 32)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(3, grid.azimuth_nodes)
    modes = [(1, 1), (3, 1), (3, 3)]
    state = ComplexF64[0.7+0.2im, 0.1-0.03im, -0.04+0.02im,
        -0.03+0.01im, 0.004+0.002im, 0.002-0.001im]
    radial = Dict((l, (m, 0)) => state[index] for (index, (l, m)) in enumerate(modes))
    slopes = Dict((l, (m, 0)) => state[3+index] for (index, (l, m)) in enumerate(modes))
    position, alpha, frequency = 8.0, 0.42, 0.4
    errors = Float64[]
    for spin in (0.0, 0.7, 0.99), potential in (:cosine, :quartic, :sextic)
        result = two_phase_observables_at_radius(radial, slopes, position, spin, alpha,
            [frequency, 0.41], [1, 2], polar, grid, -grid.azimuth_nodes, [0.0]; potential)
        reference = energy_per_tortoise_radius(state, position, frequency, spin, alpha,
            modes, polar, azimuthal, grid; potential)
        push!(errors, abs(result.energy_per_tortoise_radius-reference)/abs(reference))
        @test isapprox(result.energy_per_tortoise_radius, reference; rtol=2e-12)
    end
    println("Combined observables maximum single-cloud energy relative error: ", maximum(errors))
end

@testset "Combined observables: joint charges and currents" begin
    grid = angular_grid(16, 25)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    phase1_nodes, phase2_nodes = grid.azimuth_nodes, 2pi .* (0:26) ./ 27
    frequencies, cloud_m, alpha = [0.4, 0.41], [1, 2], 0.42
    radial = Dict((1, (1, 0)) => 0.3+0.2im, (2, (0, 1)) => -0.2+0.1im,
        (0, (2, -1)) => 0.05-0.03im, (1, (3, -2)) => 0.03+0.01im)
    slopes = Dict(pair => 0.07im*value for (pair, value) in radial)
    original = (copy(radial), copy(slopes), copy(polar[(1, 1)]))
    evaluate = (values, derivatives, position, spin; potential=:cosine) ->
        two_phase_observables_at_radius(values, derivatives, position, spin, alpha,
            frequencies, cloud_m, polar, grid, phase1_nodes, phase2_nodes; potential)

    for position in (3.0, 12.0), spin in (0.0, 0.7)
        result = evaluate(radial, slopes, position, spin)
        geometry = tortoise_geometry(position, spin)
        expected_energy_flux = 4pi * geometry.radial_factor * sum(
            sum(q .* frequencies) * imag(conj(value)*slopes[(l, q)])
            for ((l, q), value) in radial)
        expected_angular_flux = 4pi * geometry.radial_factor * sum(
            sum(q .* cloud_m) * imag(conj(value)*slopes[(l, q)])
            for ((l, q), value) in radial)
        @test result.radial_energy_flux > 0
        @test isapprox(result.radial_energy_flux, expected_energy_flux; rtol=2e-12)
        @test isapprox(result.radial_angular_momentum_flux, expected_angular_flux; rtol=2e-12)
        if spin == 0.0
            expected_momentum = 4pi * geometry.radius^2 * sum(
                sum(q .* frequencies) * sum(q .* cloud_m) * abs2(value)
                for ((l, q), value) in radial)
            @test isapprox(result.angular_momentum_per_tortoise_radius,
                expected_momentum; rtol=2e-12)
        end
    end
    baseline = evaluate(radial, slopes, 8.0, 0.7)
    inward = evaluate(radial, Dict(pair => -value for (pair, value) in slopes), 8.0, 0.7)
    @test inward.radial_energy_flux ≈ -baseline.radial_energy_flux
    @test inward.radial_angular_momentum_flux ≈ -baseline.radial_angular_momentum_flux
    @test inward.energy_per_tortoise_radius ≈ baseline.energy_per_tortoise_radius

    shifted = values -> Dict((l, q) => value*cis(-q[1]*0.37+q[2]*0.53)
        for ((l, q), value) in values)
    for potential in (:cosine, :quartic, :sextic)
        reference = evaluate(radial, slopes, 8.0, 0.7; potential)
        changed_origin = evaluate(shifted(radial), shifted(slopes), 8.0, 0.7; potential)
        for name in propertynames(reference)
            @test isapprox(getproperty(changed_origin, name), getproperty(reference, name);
                atol=1e-12, rtol=2e-12)
        end
    end

    # The exact quartic cross energy, including the polar integral, must survive
    # combining all tensor terms. Independent phases cancel the quadratic cross terms.
    first_values = Dict((1, (1, 0)) => 0.8+0.3im)
    second_values = Dict((2, (0, 1)) => 0.4-0.2im)
    energy = values -> evaluate(values, Dict(pair => 0im for pair in keys(values)),
        8.0, 0.7; potential=:quartic).energy_per_tortoise_radius
    interaction = energy(merge(first_values, second_values)) - energy(first_values) - energy(second_values)
    geometry = tortoise_geometry(8.0, 0.7)
    exact = -2pi * geometry.dr_drstar * alpha^2 * abs2(0.8+0.3im) * abs2(0.4-0.2im) *
        (9geometry.radius^2 + 0.7^2)/14
    @test isapprox(interaction, exact; rtol=2e-12)
    @test radial == original[1] && slopes == original[2] && polar[(1, 1)] == original[3]
    @test all(iszero, values(evaluate(Dict(), Dict(), 8.0, 0.7)))

    # Check differentiability of local observables, without claiming profile response.
    scaled_observables = scale -> collect(values(evaluate(
        Dict(pair => scale*value for (pair, value) in radial),
        Dict(pair => scale*value for (pair, value) in slopes), 8.0, 0.7)))
    automatic = ForwardDiff.jacobian(s -> scaled_observables(s[1]), [1.0])[:, 1]
    finite = (scaled_observables(1+1e-5) - scaled_observables(1-1e-5)) / 2e-5
    @test isapprox(automatic, finite; rtol=2e-8)
    @test_throws ArgumentError two_phase_observables_at_radius(radial, slopes, 8.0,
        0.7, alpha, frequencies, cloud_m, polar, grid, Float64[], phase2_nodes)
    bad_grid = merge(grid, (; polar_weights=grid.polar_weights[1:end-1]))
    @test_throws DimensionMismatch two_phase_observables_at_radius(radial, slopes,
        8.0, 0.7, alpha, frequencies, cloud_m, polar, bad_grid, phase1_nodes, phase2_nodes)
    println("Combined observables quartic cross energy: ", interaction, "; exact: ", exact)
end
