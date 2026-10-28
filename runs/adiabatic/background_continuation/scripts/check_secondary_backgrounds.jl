# Secondary and driven scalar responses on saved nearby primary backgrounds.
using Serialization, SHA, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

function reference_csv_row(path, primary_index)
    lines = readlines(path)
    names = Symbol.(split(first(lines), ','))
    values = only(filter(values -> parse(Int, first(values)) == primary_index,
                         split.(lines[2:end], ',')))
    return NamedTuple{Tuple(names)}(Tuple(parse.(Float64, values)))
end

function check_secondary_backgrounds(primary_checkpoint, output_directory)
    repository = normpath(joinpath(@__DIR__, "../../../.."))
    nearby = deserialize(primary_checkpoint)
    primary_path = joinpath(repository, nearby.metadata.seed_path)
    setup = load_primary_secondary_setup(primary_path)
    reference_directory = joinpath(repository, "runs/adiabatic/multimode_secondary/results")
    secondary_path = joinpath(reference_directory, "cosine_secondary.jls")
    secondary_branch = deserialize(secondary_path)
    secondary_branch.reading == "averaged" || error("require averaged secondary seed")
    branch_index = nearby.metadata.branch_index
    seed = only(filter(row -> row.primary_index == branch_index, secondary_branch.rows))
    seed.amplitude == nearby.metadata.amplitude || error("secondary seed amplitude differs")
    bytes2hex(sha256(read(primary_path))) == nearby.metadata.seed_sha256 ||
        error("primary seed changed since background solves")
    secondary_modes = setup.secondary_modes
    dissipative_modes = [(0, 0), (2, 0), (4, 0)]
    radiative_modes = [(3, 3), (5, 3), (7, 3)]
    grid = setup.grid
    polar = build_polar_values(7, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(7, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(7, polar, grid)
    tolerance = 1e-8
    far_options = (far_position=5000.0, far_step=0.25)
    reference_paths = [joinpath(reference_directory, name) for name in
        ("cosine_secondary_energy.csv", "cosine_dissipative_m0.csv", "cosine_radiative_m3.csv")]
    metadata = (;
        primary_checkpoint=relpath(abspath(primary_checkpoint), repository),
        primary_sha256=bytes2hex(sha256(read(primary_checkpoint))),
        secondary_path=relpath(secondary_path, repository),
        secondary_sha256=bytes2hex(sha256(read(secondary_path))),
        reference_sha256=[bytes2hex(sha256(read(path))) for path in reference_paths],
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        revision=strip(read(`git -C $repository rev-parse HEAD`, String)),
        branch_index, amplitude=nearby.metadata.amplitude,
        secondary_modes, dissipative_modes, radiative_modes,
        tolerance, reading=:averaged, jacobian_backend=:forwarddiff,
        far_options, n_polar=setup.primary.n_polar, n_azimuth=setup.primary.n_azimuth,
    )
    mkpath(output_directory)
    checkpoint_path = joinpath(output_directory, "secondary_backgrounds.jls")
    entries = if isfile(checkpoint_path)
        checkpoint = deserialize(checkpoint_path)
        checkpoint.metadata == metadata || error("resume metadata differs")
        checkpoint.entries
    else
        NamedTuple[]
    end
    function write_table(path, rows)
        open(path * ".tmp", "w") do stream
            println(stream, join(keys(first(rows)), ','))
            for row in rows
                println(stream, join(values(row), ','))
            end
        end
        mv(path * ".tmp", path; force=true)
    end
    compact(solution) = (; (key => value for (key, value) in pairs(solution) if key != :profiles)...)
    println("METADATA ", metadata)
    flush(stdout)
    for primary_entry in nearby.entries
        point = primary_entry.row
        if any(entry -> entry.row.label == point.label, entries)
            println("RESUME skipping ", point.label)
            continue
        end
        println("START ", point.label, " at ", now(UTC))
        flush(stdout)
        started = time()
        options = (; point.step, point.inner, point.outer, point.match)
        primary = merge(setup.primary, (; point.spin, point.alpha, options))
        trial_setup = merge(setup, (; primary, polar, azimuthal, couplings))
        solution = primary_entry.solution
        profiles = nonlinear_matching_residual(
            solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
            point.spin, point.alpha, primary.modes, couplings, polar, azimuthal, grid;
            potential=:cosine, save_profile=true, options...,
        )
        source = point.source == "saved" ?
            (; secondary=seed.solution, row=(alpha=setup.primary.alpha,)) :
            only(filter(entry -> entry.row.label == point.source, entries))
        guess = secondary_unknowns(source.secondary)
        guess[1] = complex(real(guess[1]) * point.alpha / source.row.alpha, imag(guess[1]))
        secondary = solve_secondary_point(trial_setup, profiles, guess;
            tolerance, jacobian_backend=:forwarddiff, reading=:averaged)
        println("SECONDARY ", point.label, " frequency=", secondary.frequency,
                " residual=", secondary.residual_norm, " iterations=", secondary.iterations)
        flush(stdout)
        energy = secondary_energy_and_fluxes(secondary, profiles, primary.modes,
            point.spin, point.alpha, secondary_modes, couplings, polar, azimuthal, grid)
        omega1, omega2 = real(solution.frequency), real(secondary.frequency)
        omega0, omega3 = 2omega1 - omega2, 2omega2 - omega1
        dissipative = solve_dissipative_mode(omega0, secondary, profiles, primary.modes,
            point.spin, point.alpha, secondary_modes, dissipative_modes,
            couplings, polar, azimuthal, grid)
        work0 = dissipative_source_work(dissipative, profiles, primary.modes, omega0,
            point.spin, point.alpha, secondary_modes, dissipative_modes, polar, azimuthal, grid)
        radiative = solve_radiative_mode(omega3, secondary, profiles, primary.modes,
            point.spin, point.alpha, secondary_modes, radiative_modes,
            couplings, polar, azimuthal, grid; far_options...)
        work3 = radiative_source_work(radiative, profiles, primary.modes, omega3,
            point.spin, point.alpha, secondary_modes, radiative_modes, polar, azimuthal, grid)
        row = (; point.label, point.source, point.alpha, point.spin, point.amplitude,
            primary_step=point.step, secondary_step=2point.step,
            point.inner, point.outer, point.match,
            primary_residual_change=norm(profiles.residual - solution.residual) / point.amplitude,
            primary_energy=point.energy, omega1_real=omega1,
            omega2_real=omega2, omega2_imag=imag(secondary.frequency),
            secondary_residual=secondary.residual_norm, secondary_iterations=secondary.iterations,
            secondary_energy=energy.energy, secondary_horizon_flux=energy.horizon_flux,
            secondary_outer_flux=energy.outer_flux,
            secondary_horizon_per_energy=energy.horizon_flux / energy.energy,
            secondary_balance=-(energy.horizon_flux + energy.outer_flux) /
                (2imag(secondary.frequency) * energy.energy),
            omega0, omega3, radiative_threshold_margin=omega3-point.alpha,
            m0_horizon_flux=dissipative.horizon_flux, m0_outer_flux=dissipative.outer_flux,
            m0_horizon_per_energy=dissipative.horizon_flux / energy.energy,
            m0_power_ratio=(dissipative.horizon_flux + dissipative.outer_flux) / work0,
            m0_condition=dissipative.condition_number,
            m0_relative_residual=dissipative.residual_norm / dissipative.source_residual_norm,
            m3_infinity_flux=radiative.infinity_flux, m3_horizon_flux=radiative.horizon_flux,
            m3_infinity_per_energy_squared=radiative.infinity_flux / energy.energy^2,
            m3_horizon_per_energy_squared=radiative.horizon_flux / energy.energy^2,
            m3_power_ratio=(radiative.horizon_flux + radiative.infinity_flux) / work3,
            m3_condition=radiative.condition_number,
            m3_relative_residual=radiative.residual_norm / radiative.source_residual_norm,
            elapsed_seconds=time()-started)
        push!(entries, (; row, secondary, dissipative=compact(dissipative), radiative=compact(radiative)))
        serialize(checkpoint_path * ".tmp", (; metadata, entries))
        mv(checkpoint_path * ".tmp", checkpoint_path; force=true)
        write_table(joinpath(output_directory, "secondary_backgrounds.csv"), [entry.row for entry in entries])
        println("ACCEPTED ", row)
        flush(stdout)
    end

    checks = NamedTuple[]
    function check(comparison, observable, error, limit)
        push!(checks, (; comparison, observable=String(observable), error, limit,
                       passed=isfinite(error) && error <= limit))
    end
    function compare(comparison, actual, reference, limit)
        for observable in keys(reference)
            check(comparison, observable,
                abs((getproperty(actual, observable)-getproperty(reference, observable)) /
                    getproperty(reference, observable)), limit)
        end
    end
    row_for(label) = only(filter(entry -> entry.row.label == label, entries)).row
    reference_energy, reference_m0, reference_m3 =
        [reference_csv_row(path, branch_index) for path in reference_paths]
    reference = (;
        omega2_real=real(seed.solution.frequency), omega2_imag=imag(seed.solution.frequency),
        secondary_energy=reference_energy.energy,
        secondary_horizon_per_energy=reference_energy.horizon_flux / reference_energy.energy,
        m0_horizon_per_energy=reference_m0.horizon_flux_per_energy,
        m3_infinity_per_energy_squared=reference_m3.infinity_flux_per_energy_squared,
        m3_horizon_per_energy_squared=reference_m3.horizon_flux_per_energy_squared,
    )
    compare("saved_baseline", row_for("baseline"), reference, 1e-6)
    for (label, actual_label, reference_label, limit) in (
        ("two_routes", "combined_via_spin", "combined_via_alpha", 1e-5),
        ("return_baseline", "return_baseline", "baseline", 1e-5),
        ("baseline_radial_step", "baseline_half_step", "baseline", 1e-4),
        ("combined_radial_step", "combined_half_step", "combined_via_alpha", 1e-4),
    )
        reference_row = row_for(reference_label)
        selected = NamedTuple{keys(reference)}(Tuple(getproperty(reference_row, key) for key in keys(reference)))
        compare(label, row_for(actual_label), selected, limit)
    end
    for entry in entries
        row = entry.row
        for (observable, error, limit) in (
            (:primary_residual_change, row.primary_residual_change, 1e-10),
            (:secondary_residual, row.secondary_residual, tolerance),
            (:secondary_balance, abs(row.secondary_balance-1), 1e-4),
            (:m0_power_ratio, abs(row.m0_power_ratio-1), 1e-5),
            (:m3_power_ratio, abs(row.m3_power_ratio-1), 1e-5),
            (:m0_relative_residual, row.m0_relative_residual, 1e-9),
            (:m3_relative_residual, row.m3_relative_residual, 1e-9),
            (:m0_condition, row.m0_condition, 1e8),
            (:m3_condition, row.m3_condition, 1e8),
        )
            check(row.label, observable, error, limit)
        end
        signs_correct = row.secondary_energy > 0 && row.omega2_imag > 0 &&
            row.secondary_horizon_flux < 0 && 0 < row.omega0 < row.alpha &&
            row.m0_horizon_flux > 0 && row.radiative_threshold_margin > 0 &&
            row.m3_infinity_flux > 0 && row.m3_horizon_flux < 0
        check(row.label, :channel_signs_and_thresholds, signs_correct ? 0.0 : 1.0, 0.0)
    end
    write_table(joinpath(output_directory, "secondary_checks.csv"), checks)
    for item in checks
        println("CHECK ", item)
    end
    all(item -> item.passed, checks) || error("secondary/driven check failed; inspect secondary_checks.csv")
    println("PASS: ", length(entries), " backgrounds, secondary plus m=0 and m=3, ", length(checks), " checks")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("usage: check_secondary_backgrounds.jl PRIMARY_CHECKPOINT OUTPUT_DIRECTORY")
    check_secondary_backgrounds(ARGS...)
end
