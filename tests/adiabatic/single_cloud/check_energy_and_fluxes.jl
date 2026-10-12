using Test, Serialization
include(joinpath(@__DIR__, "../../../src/adiabatic/load.jl"))

@testset "Energy and fluxes reproduce a saved solution" begin
    checkpoint_path = joinpath(
        @__DIR__,
        "../../../runs/adiabatic/amplitude_evolution/results/",
        "cosine_derivatives_part2.jls",
    )
    saved = deserialize(checkpoint_path)
    settings = saved.metadata

    # Branch point 51 is the previously checked A0 ≈ 2.565 solution.
    reference = only(filter(row -> row.branch_index == 51, saved.rows))
    entry = only(filter(row -> row.branch_index == 51, saved.solutions))
    amplitude = entry.amplitude
    solution = entry.solution

    @test amplitude == reference.amplitude
    @test reference.residual <= settings.tolerance

    maximum_l = maximum(first, settings.modes)
    grid = angular_grid(settings.n_polar, settings.n_azimuth)
    polar_values = build_polar_values(maximum_l, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(
        maximum_l, polar_values, grid,
    )

    packed_state = pack_complex(vcat(
        solution.frequency,
        solution.horizon_amplitudes ./ amplitude,
        solution.outer_amplitudes[2:end] ./ amplitude,
    ))

    observables = (state, trial_amplitude) -> energy_and_fluxes(
        state, trial_amplitude,
        settings.spin, settings.alpha, settings.modes, couplings,
        polar_values, azimuthal_values, grid;
        settings.options...,
    )

    computed = observables(packed_state, amplitude)

    expected = [
        reference.energy,
        reference.horizon_flux,
        reference.outer_flux,
    ]

    for (name, actual, target) in zip(
        ("energy", "horizon_flux", "outer_flux"), computed, expected,
    )
        println(name, ": computed=", actual, " saved=", target)
        @test isapprox(actual, target; rtol=1e-10, atol=0)
    end

        residual = (state, trial_amplitude) -> pack_complex(
        scaled_nonlinear_matching_residual(
            unpack_complex(state), trial_amplitude,
            settings.spin, settings.alpha, settings.modes, couplings,
            polar_values, azimuthal_values, grid;
            settings.options...,
        ),
    )

    derivatives = implicit_observable_derivatives(
        residual, observables, packed_state, amplitude,
    )
    slopes = derivatives.observable_derivatives

    @test all(isfinite, slopes)
    @test isapprox(
        slopes[1], reference.energy_derivative; rtol=1e-9, atol=0,
    )
    println("Implicit slopes [E′, F_H′, F_out′]: ", slopes)
    flush(stdout)

    function corrected_observables(trial_amplitude)
        # Keep the baseline normalized coefficients as the initial guess.
        # Do not use the implicit derivative being tested.
        amplitude_ratio = trial_amplitude / amplitude

        corrected_solution = solve_nonlinear_mode(
            solution.frequency, trial_amplitude,
            settings.spin, settings.alpha, settings.modes, couplings,
            polar_values, azimuthal_values, grid;
            horizon_guess =
                amplitude_ratio .* solution.horizon_amplitudes,
            outer_tail_guess =
                amplitude_ratio .* solution.outer_amplitudes[2:end],
            tolerance=settings.tolerance,
            max_iterations=12,
            jacobian_backend=:forwarddiff,
            settings.options...,
        )

        @test corrected_solution.scaled_residual_norm <= settings.tolerance

        corrected_state = pack_complex(vcat(
            corrected_solution.frequency,
            corrected_solution.horizon_amplitudes ./ trial_amplitude,
            corrected_solution.outer_amplitudes[2:end] ./ trial_amplitude,
        ))

        return observables(corrected_state, trial_amplitude)
    end

    error_history = Vector{Float64}[]

    for fraction in (0.01, 0.005, 0.0025)
        spacing = fraction * amplitude
        above = corrected_observables(amplitude + spacing)
        below = corrected_observables(amplitude - spacing)

        finite_difference = (above .- below) ./ (2spacing)
        relative_errors = abs.((finite_difference .- slopes) ./ slopes)
        push!(error_history, relative_errors)

        println("h=", spacing,
                "  finite-difference slopes=", finite_difference,
                "  relative errors=", relative_errors)
        flush(stdout)
    end

    # Check each observable separately at this branch point.
    for index in 1:3
        @test error_history[end][index] < 1e-4
        @test error_history[end][index] < error_history[1][index]
    end
end
