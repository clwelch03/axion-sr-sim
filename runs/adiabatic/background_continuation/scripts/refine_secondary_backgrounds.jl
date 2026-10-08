# Resolve the two m=3 horizon-flux refinement failures without changing their
# original results or tolerance. Re-solve primary, secondary, and driven fields
# at primary step 0.025, comparing against the saved step-0.05 calculations.
using Serialization, SHA, Dates
include(joinpath(@__DIR__, "check_secondary_backgrounds.jl"))

function refine_secondary_backgrounds()
    repository = normpath(joinpath(@__DIR__, "../../../.."))
    primary_path = joinpath(@__DIR__, "../results/amplitude2052/nearby_backgrounds.jls")
    response_path = joinpath(@__DIR__, "../results/amplitude2052/secondary/secondary_backgrounds.jls")
    nearby, responses = deserialize(primary_path), deserialize(response_path)
    setup = load_secondary_setup(joinpath(repository, nearby.metadata.seed_path))
    primary_modes, secondary_modes = setup.primary.modes, setup.secondary_modes
    dissipative_modes, radiative_modes = responses.metadata.dissipative_modes, responses.metadata.radiative_modes
    grid = setup.grid
    polar = build_polar_values(7, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(7, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(7, polar, grid)
    metadata = (; primary_sha256=bytes2hex(sha256(read(primary_path))),
        response_sha256=bytes2hex(sha256(read(response_path))),
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        revision=strip(read(`git -C $repository rev-parse HEAD`, String)),
        primary_step=0.025, primary_tolerance=1e-11, secondary_tolerance=1e-8,
        parent_metadata=responses.metadata)
    output = joinpath(@__DIR__, "../results/amplitude2052/secondary/quarter_step")
    mkpath(output)
    checkpoint_path = joinpath(output, "refined_responses.jls")
    entries = if isfile(checkpoint_path)
        saved = deserialize(checkpoint_path)
        saved.metadata == metadata || error("refinement resume metadata differs")
        saved.entries
    else
        NamedTuple[]
    end
    function write_table(path, rows)
        open(path * ".tmp", "w") do io
            println(io, join(keys(first(rows)), ','))
            for row in rows
                println(io, join(values(row), ','))
            end
        end
        mv(path * ".tmp", path; force=true)
    end
    compact(solution) = (; (key => value for (key, value) in pairs(solution) if key != :profiles)...)
    println("METADATA ", metadata)
    flush(stdout)
    for label in ("baseline_half_step", "combined_half_step")
        any(entry -> entry.row.label == label, entries) && continue
        source = only(filter(entry -> entry.row.label == label, nearby.entries))
        previous = only(filter(entry -> entry.row.label == label, responses.entries))
        point = source.row
        println("START ", label, " at ", now(UTC))
        flush(stdout)
        started = time()
        options = (; step=metadata.primary_step, point.inner, point.outer, point.match)
        primary = solve_nonlinear_mode(source.solution.frequency, point.amplitude,
            point.spin, point.alpha, primary_modes, couplings, polar, azimuthal, grid;
            horizon_guess=source.solution.horizon_amplitudes,
            outer_tail_guess=source.solution.outer_amplitudes[2:end], potential=:cosine,
            jacobian_backend=:forwarddiff, tolerance=metadata.primary_tolerance, options...)
        profiles = nonlinear_matching_residual(primary.frequency, primary.horizon_amplitudes,
            primary.outer_amplitudes, point.spin, point.alpha, primary_modes,
            couplings, polar, azimuthal, grid; potential=:cosine, save_profile=true, options...)
        secondary = solve_secondary_mode(previous.secondary.frequency, profiles, primary_modes,
            point.spin, point.alpha, secondary_modes, couplings, polar, azimuthal, grid;
            horizon_guess=previous.secondary.horizon_amplitudes,
            outer_tail_guess=previous.secondary.outer_amplitudes[2:end],
            jacobian_backend=:forwarddiff, reading=:averaged, tolerance=metadata.secondary_tolerance)
        energy = secondary_energy_and_fluxes(secondary, profiles, primary_modes,
            point.spin, point.alpha, secondary_modes, couplings, polar, azimuthal, grid)
        omega0 = 2real(primary.frequency)-real(secondary.frequency)
        omega3 = 2real(secondary.frequency)-real(primary.frequency)
        dissipative = solve_dissipative_mode(omega0, secondary, profiles, primary_modes,
            point.spin, point.alpha, secondary_modes, dissipative_modes, couplings, polar, azimuthal, grid)
        work0 = dissipative_source_work(dissipative, profiles, primary_modes, omega0,
            point.spin, point.alpha, secondary_modes, dissipative_modes, polar, azimuthal, grid)
        radiative = solve_radiative_mode(omega3, secondary, profiles, primary_modes,
            point.spin, point.alpha, secondary_modes, radiative_modes, couplings, polar, azimuthal, grid;
            responses.metadata.far_options...)
        work3 = radiative_source_work(radiative, profiles, primary_modes, omega3,
            point.spin, point.alpha, secondary_modes, radiative_modes, polar, azimuthal, grid)
        row = (; label, point.alpha, point.spin, point.amplitude, primary_step=metadata.primary_step,
            primary_residual=primary.scaled_residual_norm, secondary_residual=secondary.residual_norm,
            omega2_real=real(secondary.frequency), omega2_imag=imag(secondary.frequency),
            secondary_energy=energy.energy, secondary_horizon_per_energy=energy.horizon_flux/energy.energy,
            secondary_balance=-(energy.horizon_flux+energy.outer_flux)/(2imag(secondary.frequency)*energy.energy),
            m0_horizon_per_energy=dissipative.horizon_flux/energy.energy,
            m0_power_ratio=(dissipative.horizon_flux+dissipative.outer_flux)/work0,
            m3_infinity_per_energy_squared=radiative.infinity_flux/energy.energy^2,
            m3_horizon_per_energy_squared=radiative.horizon_flux/energy.energy^2,
            m3_power_ratio=(radiative.horizon_flux+radiative.infinity_flux)/work3,
            m0_relative_residual=dissipative.residual_norm/dissipative.source_residual_norm,
            m3_relative_residual=radiative.residual_norm/radiative.source_residual_norm,
            elapsed_seconds=time()-started)
        push!(entries, (; row, primary, secondary, dissipative=compact(dissipative), radiative=compact(radiative)))
        serialize(checkpoint_path * ".tmp", (; metadata, entries))
        mv(checkpoint_path * ".tmp", checkpoint_path; force=true)
        write_table(joinpath(output, "refined_responses.csv"), [entry.row for entry in entries])
        println("ACCEPTED ", row)
        flush(stdout)
    end
    checks = NamedTuple[]
    function check(label, observable, error, limit)
        push!(checks, (; label, observable=String(observable), error, limit,
                       passed=isfinite(error) && error <= limit))
    end
    for entry in entries
        row = entry.row
        previous = only(filter(entry -> entry.row.label == row.label, responses.entries)).row
        coarse_label = row.label == "baseline_half_step" ? "baseline" : "combined_via_alpha"
        coarse = only(filter(entry -> entry.row.label == coarse_label, responses.entries)).row
        for key in (:omega2_real, :omega2_imag, :secondary_energy, :secondary_horizon_per_energy,
                    :m0_horizon_per_energy, :m3_infinity_per_energy_squared, :m3_horizon_per_energy_squared)
            value, prior = getproperty(row, key), getproperty(previous, key)
            change = abs((value-prior)/prior)
            check(row.label, key, change, 1e-4)
            if key == :m3_horizon_per_energy_squared
                reduction = abs((getproperty(coarse,key)-prior)/(prior-value))
                println("REFINEMENT ", row.label, " m3 horizon relative change=", change,
                        " reduction factor=", reduction, " (16 expected for leading RK4 error)")
            end
        end
        for (key, error, limit) in (
            (:primary_residual, row.primary_residual, metadata.primary_tolerance),
            (:secondary_residual, row.secondary_residual, metadata.secondary_tolerance),
            (:secondary_balance, abs(row.secondary_balance-1), 1e-4),
            (:m0_power_ratio, abs(row.m0_power_ratio-1), 1e-5),
            (:m3_power_ratio, abs(row.m3_power_ratio-1), 1e-5),
            (:m0_relative_residual, row.m0_relative_residual, 1e-9),
            (:m3_relative_residual, row.m3_relative_residual, 1e-9),
        )
            check(row.label, key, error, limit)
        end
    end
    write_table(joinpath(output, "refinement_checks.csv"), checks)
    all(check -> check.passed, checks) || error("quarter-step refinement failed; inspect refinement_checks.csv")
    println("PASS: ", length(entries), " refined backgrounds, ", length(checks), " checks")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    refine_secondary_backgrounds()
end
