using Serialization, LinearAlgebra

"""
    single_cloud_continuation_entry_adapter(point_id, entry; tolerance)

Convert one accepted single-cloud entry, with fields `amplitude` and `solution`,
to a ContinuationPoint. Supply the tolerance used for that solve; the entry
does not store it. Use `nothing` if the original tolerance is unknown.

Copy the boundary coefficients and residual so the new record owns its arrays.
This conversion does not solve, check convergence, or reconstruct radial samples.
"""
function single_cloud_continuation_entry_adapter(point_id, legacy_continuation_record; tolerance)
    (; amplitude, solution) = legacy_continuation_record

    # A single nonlinear cloud has one base frequency, including when its
    # profile contains several angular modes and higher harmonics.
    profile = TwoPhaseProfile(
        [solution.frequency],
        copy(solution.horizon_amplitudes),
        copy(solution.outer_amplitudes),
    )

    # The imposed parameter is the first physical outer coefficient, A0.
    parameter_values = (; amplitude)

    # Newton tests norm(residual / amplitude) against tolerance.
    # Preserve the stored norms rather than recomputing them during conversion.
    diagnostics = (;
        residual=copy(solution.residual),
        residual_norm=solution.residual_norm,
        scaled_residual_norm=solution.scaled_residual_norm,
        residual_scale=amplitude,
        tolerance,
        iterations=solution.iterations,
    )

    # This model has no fixed cloud background, and these entries contain
    # boundary data only, with no saved radial integration samples.
    backgrounds = ()
    radial_result = nothing

    return ContinuationPoint(
        point_id, parameter_values, profile, diagnostics, backgrounds, radial_result,
    )
end

"""
    single_cloud_continuation_checkpoint_adapter(legacy_checkpoint; checkpoint_path)

Convert an old single-cloud checkpoint containing `continuation` and `rows`
to a ContinuationBranch. Preserve the other saved fields as branch settings
and each saved table row as `point.diagnostics.observables`.

The file records a branch-level tolerance, not a tolerance for every past solve.
Keep that setting, but leave each point's original tolerance as `nothing`.
This conversion does not solve or recheck convergence.
"""
function single_cloud_continuation_checkpoint_adapter(legacy_checkpoint; checkpoint_path)
    legacy_continuation_records = legacy_checkpoint.continuation
    legacy_diagnostic_rows = legacy_checkpoint.rows
    length(legacy_continuation_records) == length(legacy_diagnostic_rows) ||
        throw(ArgumentError("Require one diagnostic row per saved continuation point."))

    # Keep the saved physical, grid, Newton, and continuation settings.
    # The two histories are represented by the points below.
    branch_settings = (;
        model=:single_cloud,
        (field_name => deepcopy(field_value)
         for (field_name, field_value) in pairs(legacy_checkpoint)
         if field_name ∉ (:continuation, :rows))...,
    )

    branch_points = map(eachindex(legacy_continuation_records)) do point_id
        legacy_continuation_record = legacy_continuation_records[point_id]
        legacy_diagnostic_row = legacy_diagnostic_rows[point_id]
        legacy_diagnostic_row.amplitude == legacy_continuation_record.amplitude ||
            throw(ArgumentError("Diagnostic row $point_id has a different amplitude."))

        # Original positions become stable IDs when importing the old file.
        converted_point = single_cloud_continuation_entry_adapter(
            point_id, legacy_continuation_record; tolerance=nothing,
        )
        point_diagnostics = (;
            converted_point.diagnostics...,
            observables=deepcopy(legacy_diagnostic_row),
        )

        return ContinuationPoint(
            converted_point.point_id,
            converted_point.parameter_values,
            converted_point.profile,
            point_diagnostics,
            converted_point.backgrounds,
            converted_point.radial_result,
        )
    end

    provenance = (;
        checkpoint_path=abspath(checkpoint_path),
        source_format=:legacy_single_cloud,
    )
    return ContinuationBranch(branch_settings, branch_points, provenance)
end

"""
    load_branch_setup(checkpoint_path; potential)

Read a branch checkpoint and build the angular grid and couplings for its
settings. Return the checkpoint (`saved`), its physical and radial settings,
the angular data, and the requested `potential`.
"""
function load_branch_setup(checkpoint_path; potential)
    potential = Symbol(potential)
    potential in (:cosine, :quartic, :sextic) ||
        throw(ArgumentError("potential must be cosine, quartic, or sextic"))

    saved = deserialize(checkpoint_path)
    maximum_l = maximum(first, saved.modes)
    grid = angular_grid(saved.n_polar, saved.n_azimuth)
    polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_l, polar, grid)

    return (;
        checkpoint_path, saved, potential,
        spin=saved.spin, alpha=saved.alpha, modes=saved.modes, options=saved.options,
        n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
        grid, polar, azimuthal, couplings,
    )
end

# Finite-radius outward current, including the conjugate-mode normalization.
function outward_flux(setup, state, position, frequency)
    count = length(setup.modes)
    geometry = tortoise_geometry(position, setup.spin)

    return 4π * geometry.radial_factor * sum(
        real(im * m * frequency * state[index] * conj(state[count + index]))
        for (index, (_, m)) in enumerate(setup.modes)
    )
end


"""
    branch_diagnostics(setup, amplitude, solution)

Reintegrate an accepted solution and return its energy, signed fluxes,
frequency, peak field, and convergence data as one table row.
"""
function branch_diagnostics(setup, amplitude, solution)
    (; spin, alpha, modes, couplings, polar, azimuthal, grid, options, potential) = setup

    # Retain profiles for the energy integral.
    profiles = nonlinear_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        potential, save_profile=true, options...,
    )
    branches = (profiles.from_horizon, profiles.from_outer)

    energy = sum(
        integrate_energy_branch(
            branch, solution.frequency, spin, alpha, modes, polar, azimuthal, grid;
            potential,
        )
        for branch in branches
    )

    # Loss into the horizon has the opposite sign to the outward current there.
    horizon_flux = -outward_flux(
        setup, profiles.from_horizon.states[:, 1], options.inner, solution.frequency,
    )
    outer_flux = outward_flux(
        setup, profiles.from_outer.states[:, 1], options.outer, solution.frequency,
    )
    total_flux = horizon_flux + outer_flux

    # Estimate the peak field using every tenth radial sample.
    peak_field = 0.0
    for branch in branches, index in 1:10:length(branch.positions)
        radial_values = Dict(
            mode => branch.states[k, index] for (k, mode) in enumerate(modes)
        )
        field = reconstruct_field(radial_values, polar, azimuthal, grid)
        peak_field = max(peak_field, maximum(abs, field))
    end

    return (;
        amplitude, energy, normalized_energy=energy / amplitude^2,
        omega_real=real(solution.frequency), omega_imag=imag(solution.frequency),
        horizon_flux, outer_flux, total_flux,
        balance_ratio=-total_flux / (2 * imag(solution.frequency) * energy),
        peak_field, residual=solution.scaled_residual_norm,
        iterations=solution.iterations,
    )
end

"""
    save_branch(output_prefix, setup, settings, continuation, rows)

Write the full solutions (`.jls`, for restarting) and the table (`.csv`),
each replaced only after its temporary file is complete.
"""
function save_branch(output_prefix, setup, settings, continuation, rows)
    data = (;
        setup.spin, setup.alpha, setup.modes, setup.options,
        setup.n_polar, setup.n_azimuth, potential=String(setup.potential),
        settings.tolerance, predictor=String(settings.predictor),
        settings.history_spacing, jacobian_backend=String(settings.jacobian_backend),
        settings.maximum_amplitude_step, settings.maximum_fractional_step,
        settings.newton_max_iterations, settings.difference_step,
        continuation, rows,
    )
    serialize(output_prefix * ".jls.tmp", data)
    mv(output_prefix * ".jls.tmp", output_prefix * ".jls"; force=true)

    open(output_prefix * ".csv.tmp", "w") do stream
        println(stream, join(string.(keys(first(rows))), ','))
        for row in rows
            println(stream, join(values(row), ','))
        end
    end
    mv(output_prefix * ".csv.tmp", output_prefix * ".csv"; force=true)
end

"""
    start_branch(seed_path, output_prefix, potential, settings)

Choose where a continuation starts. Resume `output_prefix.jls` when it exists.
Otherwise start from the seed: a seed with the requested potential supplies
its whole history; a seed with another potential supplies only its first
solution, which is re-solved as a guess for the requested potential.
"""
function start_branch(seed_path, output_prefix, potential, settings)
    resuming = isfile(output_prefix * ".jls")
    setup = load_branch_setup(resuming ? output_prefix * ".jls" : seed_path; potential)
    saved = setup.saved

    if Symbol(saved.potential) == setup.potential
        return (; setup, resuming,
            continuation=copy(saved.continuation), rows=copy(saved.rows))
    end
    resuming && error("checkpoint potential does not match the requested potential")

    (; spin, alpha, modes, couplings, polar, azimuthal, grid, options) = setup
    seed_entry = first(saved.continuation)
    seed_solution = seed_entry.solution
    solution = solve_nonlinear_mode(
        seed_solution.frequency, seed_entry.amplitude,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        horizon_guess=seed_solution.horizon_amplitudes,
        outer_tail_guess=seed_solution.outer_amplitudes[2:end],
        potential=setup.potential,
        settings.tolerance, max_iterations=settings.newton_max_iterations,
        settings.difference_step, settings.jacobian_backend, options...,
    )
    row = branch_diagnostics(setup, seed_entry.amplitude, solution)
    println("SEED ", row)
    flush(stdout)

    return (; setup, resuming,
        continuation=[(amplitude=seed_entry.amplitude, solution)], rows=[row])
end

"""
    continue_branch(seed_path, output_prefix; potential, target_amplitude, kwargs...)

Continue a stationary branch in increasing amplitude until `target_amplitude`
or `runtime_seconds`, saving after every accepted point.

Before each attempt the step is limited to
`min(maximum_amplitude_step, maximum_fractional_step * A0)`. A failed
prediction or correction halves the step; cheap corrections grow it.
"""
function continue_branch(
    seed_path, output_prefix;
    potential, target_amplitude,
    runtime_seconds=1000.0,
    predictor=:secant, history_spacing=0.025,
    jacobian_backend=:forwarddiff, difference_step=1e-8,
    newton_max_iterations=12, tolerance=1e-9,
    maximum_amplitude_step=0.25, maximum_fractional_step=0.05,
    initial_step=nothing,
)
    predictor in (:secant, :cubic) ||
        throw(ArgumentError("predictor must be secant or cubic"))
    for (name, value) in (
        (:target_amplitude, target_amplitude), (:runtime_seconds, runtime_seconds),
        (:history_spacing, history_spacing), (:difference_step, difference_step),
        (:maximum_amplitude_step, maximum_amplitude_step),
        (:maximum_fractional_step, maximum_fractional_step),
    )
        value > 0 || throw(ArgumentError("require positive $name"))
    end
    newton_max_iterations > 0 ||
        throw(ArgumentError("require a positive correction limit"))

    settings = (;
        tolerance, predictor, history_spacing, jacobian_backend, difference_step,
        newton_max_iterations, maximum_amplitude_step, maximum_fractional_step,
    )
    started = time()

    (; setup, resuming, continuation, rows) =
        start_branch(seed_path, output_prefix, potential, settings)
    (; spin, alpha, modes, couplings, polar, azimuthal, grid, options) = setup
    println("potential=", setup.potential, " resuming=", resuming,
        " start=", setup.checkpoint_path, " output=", output_prefix)
    println("SETTINGS ", settings, " target_amplitude=", target_amplitude)
    flush(stdout)

    # Respect small accepted steps when restarting a difficult branch.
    amplitude = last(continuation).amplitude
    amplitude_step = 0.1 * amplitude
    if length(continuation) >= 2
        previous_step = amplitude - continuation[end - 1].amplitude
        amplitude_step = min(1.5 * previous_step, amplitude_step)
    end
    if !isnothing(initial_step)
        initial_step > 0 && isfinite(initial_step) ||
            throw(ArgumentError("require positive finite initial step"))
        amplitude_step = initial_step
    end

    save_branch(output_prefix, setup, settings, continuation, rows)

    progress = function(iteration, residual_norm)
        println("  Newton iteration=$iteration scaled_residual=$residual_norm")
        flush(stdout)
    end

    while time() - started < runtime_seconds && amplitude < target_amplitude
        # Enforce both limits for initial, restarted, grown, and retried steps.
        amplitude_step = min(amplitude_step, maximum_amplitude_step,
                             maximum_fractional_step * amplitude)
        # Predict the next point without advancing beyond the requested target.
        target = min(target_amplitude, amplitude + amplitude_step)
        prediction = if predictor === :cubic
            predict_nonlinear_mode_cubic(continuation, target; history_spacing)
        else
            predict_nonlinear_mode(continuation, target)
        end

        valid_frequency = isfinite(prediction.frequency) &&
            0 < real(prediction.frequency) < alpha
        valid_amplitudes = all(isfinite, prediction.horizon_amplitudes) &&
            all(isfinite, prediction.outer_amplitudes)

        if !valid_frequency || !valid_amplitudes
            amplitude_step /= 2
            amplitude_step < 1e-5 * amplitude &&
                error("No admissible continuation prediction")
            continue
        end

        println("Starting A0=$target (elapsed $(round(time() - started)) s)")
        flush(stdout)

        # Correct the guess. On failure, retry from the last accepted point.
        solution = try
            solve_nonlinear_mode(
                prediction.frequency, target,
                spin, alpha, modes, couplings, polar, azimuthal, grid;
                horizon_guess=prediction.horizon_amplitudes,
                outer_tail_guess=prediction.outer_amplitudes[2:end],
                potential=setup.potential,
                tolerance, max_iterations=newton_max_iterations,
                difference_step, jacobian_backend,
                progress_callback=progress, options...,
            )
        catch exception
            exception isa ErrorException || rethrow()
            println("FAILED A0=$target: ", sprint(showerror, exception))
            flush(stdout)

            amplitude_step /= 2
            amplitude_step < 1e-5 * amplitude &&
                error("Continuation step exhausted; last accepted A0=$amplitude")
            continue
        end

        # Measure and save the accepted point before attempting another.
        row = branch_diagnostics(setup, target, solution)
        push!(continuation, (; amplitude=target, solution))
        push!(rows, row)
        save_branch(output_prefix, setup, settings, continuation, rows)

        println("SAVED ", row, " elapsed=", round(time() - started))
        flush(stdout)

        # Cheap corrections allow a larger next step; expensive ones reduce it.
        accepted_step = target - amplitude
        amplitude = target

        growth = if solution.iterations <= 4
            1.6
        elseif solution.iterations <= 8
            1.3
        elseif solution.iterations <= 16
            1.0
        else
            0.75
        end

        amplitude_step = growth * accepted_step
    end

    println("Finished at A0=", amplitude, " with ", length(rows), " saved points.")
    return (; setup, continuation, rows)
end
