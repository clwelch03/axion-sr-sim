# Stationary continuation in alpha and spin at fixed outer amplitude.
# Reuse the existing solver; do not evolve the BH or either cloud population.
using Serialization, SHA, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/single_cloud/load.jl"))

function check_nearby_backgrounds(output_directory)
    repository = normpath(joinpath(@__DIR__, "../../../.."))
    seed_path = joinpath(repository, "runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
    setup = load_branch_setup(seed_path; potential=:cosine)
    seed_entry = first(setup.saved.continuation)
    amplitude = seed_entry.amplitude
    tolerance = 1e-11
    metadata = (;
        seed_path=relpath(seed_path, repository), seed_sha256=bytes2hex(sha256(read(seed_path))),
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        revision=strip(read(`git -C $repository rev-parse HEAD`, String)),
        amplitude, tolerance, jacobian_backend=:forwarddiff, max_iterations=12,
        setup.modes, setup.n_polar, setup.n_azimuth, potential=:cosine,
        base_options=setup.options, base_alpha=setup.alpha, base_spin=setup.spin,
        alpha_spacing=0.001, spin_spacing=0.001,
    )
    mkpath(output_directory)
    checkpoint_path = joinpath(output_directory, "nearby_backgrounds.jls")
    entries = if isfile(checkpoint_path)
        checkpoint = deserialize(checkpoint_path)
        checkpoint.metadata == metadata || error("resume metadata differs from this run")
        checkpoint.entries
    else
        NamedTuple[]
    end
    println("METADATA ", metadata)
    flush(stdout)

    # Two routes to the joint offset, a return to the original background,
    # and actual re-solves at half the radial spacing.
    cases = [
        (label="baseline", source="saved", alpha=0.42, spin=0.99, step=0.1),
        (label="alpha_minus", source="baseline", alpha=0.419, spin=0.99, step=0.1),
        (label="alpha_plus", source="baseline", alpha=0.421, spin=0.99, step=0.1),
        (label="spin_minus", source="baseline", alpha=0.42, spin=0.989, step=0.1),
        (label="spin_plus", source="baseline", alpha=0.42, spin=0.991, step=0.1),
        (label="combined_via_alpha", source="alpha_minus", alpha=0.419, spin=0.989, step=0.1),
        (label="combined_via_spin", source="spin_minus", alpha=0.419, spin=0.989, step=0.1),
        (label="return_baseline", source="combined_via_alpha", alpha=0.42, spin=0.99, step=0.1),
        (label="baseline_half_step", source="baseline", alpha=0.42, spin=0.99, step=0.05),
        (label="combined_half_step", source="combined_via_alpha", alpha=0.419, spin=0.989, step=0.05),
    ]
    setup.alpha == 0.42 && setup.spin == 0.99 && setup.options.step == 0.1 ||
        error("unexpected reference background or radial step")

    function save_results()
        serialize(checkpoint_path * ".tmp", (; metadata, entries))
        mv(checkpoint_path * ".tmp", checkpoint_path; force=true)
        table_path = joinpath(output_directory, "nearby_backgrounds.csv")
        open(table_path * ".tmp", "w") do stream
            println(stream, join(keys(first(entries).row), ','))
            for entry in entries
                println(stream, join(values(entry.row), ','))
            end
        end
        mv(table_path * ".tmp", table_path; force=true)
    end

    for case in cases
        if any(entry -> entry.row.label == case.label, entries)
            println("RESUME skipping ", case.label)
            continue
        end
        source = case.source == "saved" ?
            (; solution=seed_entry.solution, row=(alpha=setup.alpha,)) :
            only(filter(entry -> entry.row.label == case.source, entries))
        trial_setup = merge(setup, (; case.alpha, case.spin,
            options=merge(setup.options, (; case.step))))
        # Preserve omega_R / alpha only as a predictor. Newton corrects omega
        # and every free boundary coefficient on the target Kerr geometry.
        frequency_guess = complex(real(source.solution.frequency) * case.alpha / source.row.alpha,
                                  imag(source.solution.frequency))
        println("START ", case, " at ", now(UTC))
        flush(stdout)
        started = time()
        solution = solve_nonlinear_mode(
            frequency_guess, amplitude, case.spin, case.alpha,
            setup.modes, setup.couplings, setup.polar, setup.azimuthal, setup.grid;
            horizon_guess=source.solution.horizon_amplitudes,
            outer_tail_guess=source.solution.outer_amplitudes[2:end],
            potential=:cosine, tolerance, max_iterations=12, jacobian_backend=:forwarddiff,
            progress_callback=(iteration, residual) -> begin
                println("NEWTON ", case.label, " iteration=", iteration, " residual=", residual)
                flush(stdout)
            end,
            trial_setup.options...,
        )
        diagnostics = branch_diagnostics(trial_setup, amplitude, solution)
        # The retained primary harmonics all have m_n/omega_n = 1/omega_R.
        horizon_angular_velocity = case.spin / (2 * (1 + sqrt(1 - case.spin^2)))
        row = (; case.label, source=case.source, case.alpha, case.spin, case.step,
            setup.options.inner, setup.options.outer, setup.options.match,
            diagnostics...,
            horizon_torque=diagnostics.horizon_flux / diagnostics.omega_real,
            superradiant_margin=horizon_angular_velocity - diagnostics.omega_real,
            elapsed_seconds=time() - started)
        push!(entries, (; row, solution))
        save_results()
        println("ACCEPTED ", row)
        flush(stdout)
    end

    row_for(label) = only(filter(entry -> entry.row.label == label, entries)).row
    checks = NamedTuple[]
    function compare_rows(label, actual, reference; relative_tolerance)
        for observable in (:omega_real, :omega_imag, :energy, :horizon_flux)
            actual_value, reference_value = getproperty(actual, observable), getproperty(reference, observable)
            relative_error = abs((actual_value - reference_value) / reference_value)
            push!(checks, (; comparison=label, observable=String(observable), relative_error,
                           relative_tolerance, passed=relative_error <= relative_tolerance))
        end
    end
    compare_rows("saved_baseline", row_for("baseline"), first(setup.saved.rows); relative_tolerance=1e-7)
    compare_rows("two_routes", row_for("combined_via_spin"), row_for("combined_via_alpha"); relative_tolerance=1e-5)
    compare_rows("return_baseline", row_for("return_baseline"), row_for("baseline"); relative_tolerance=1e-5)
    compare_rows("baseline_radial_step", row_for("baseline_half_step"), row_for("baseline"); relative_tolerance=1e-4)
    compare_rows("combined_radial_step", row_for("combined_half_step"), row_for("combined_via_alpha"); relative_tolerance=1e-4)
    open(joinpath(output_directory, "checks.csv"), "w") do stream
        println(stream, join(keys(first(checks)), ','))
        for check in checks
            println(stream, join(values(check), ','))
            println("CHECK ", check)
        end
    end
    all(entry -> entry.row.residual <= tolerance, entries) || error("unconverged matching residual")
    all(entry -> entry.row.horizon_flux < 0 && entry.row.superradiant_margin > 0, entries) ||
        error("unexpected loss of superradiance in this local neighborhood")
    all(check -> check.passed, checks) || error("stationary comparison failed; inspect checks.csv")
    println("PASS: ", length(entries), " stationary solves and ", length(checks), " observable comparisons")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("usage: check_nearby_backgrounds.jl [OUTPUT_DIRECTORY]")
    check_nearby_backgrounds(isempty(ARGS) ? joinpath(@__DIR__, "../results") : only(ARGS))
end
