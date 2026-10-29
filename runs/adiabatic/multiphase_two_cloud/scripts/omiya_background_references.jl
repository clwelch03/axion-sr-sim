# Continue the multimode Omiya reference to a new Kerr spin, then match E1/E2.
# Usage: julia --project=. omiya_background_references.jl INPUT_DIRECTORY OUTPUT_DIRECTORY
include("omiya_energy_references.jl")

function correct_reference_primary(setup, previous, amplitude, prefix)
    p, old = setup.primary, previous.solution
    ratio = amplitude / previous.amplitude
    progress = (i, residual) -> begin
        println("PRIMARY spin=", p.spin, " amplitude=", amplitude,
                " Newton=", i, " residual=", residual); flush(stdout)
    end
    solution = solve_nonlinear_mode(old.frequency, amplitude, p.spin, p.alpha,
        p.modes, setup.couplings, setup.polar, setup.azimuthal, setup.grid;
        horizon_guess=ratio .* old.horizon_amplitudes,
        outer_tail_guess=ratio .* old.outer_amplitudes[2:end],
        potential=Symbol(p.potential), p.options..., tolerance=1e-9,
        jacobian_backend=:forwarddiff, max_iterations=20, progress_callback=progress)
    profiles = reconstruct_primary(setup, solution)
    energy = primary_energy(setup, solution, profiles)
    @assert norm(profiles.residual) / amplitude < 1e-8
    @assert isfinite(energy) && energy > 0
    @assert all(isfinite, profiles.from_horizon.states) && all(isfinite, profiles.from_outer.states)
    current = (; amplitude, solution, profiles, energy, previous.seed_indices)
    serialize(prefix * ".jls", (; current, spin=p.spin, alpha=p.alpha))
    println("PRIMARY ENERGY spin=", p.spin, " energy=", energy); flush(stdout)
    current
end

function match_reference_energy(setup, current, target, prefix)
    previous_x = previous_error = nothing
    derivative = 2.0
    for iteration in 1:8
        abs(current.energy / target - 1) < 1e-6 && return current
        x, energy_error = log(current.amplitude), log(current.energy / target)
        if previous_x !== nothing
            derivative = (energy_error - previous_error) / (x - previous_x)
            isfinite(derivative) && derivative > 0 || error("Nonmonotone local energy response")
        end
        next_x = x + clamp(-energy_error / derivative, -0.1, 0.1)
        previous_x, previous_error = x, energy_error
        current = correct_reference_primary(setup, current, exp(next_x), prefix * "-step$iteration")
    end
    abs(current.energy / target - 1) < 1e-6 || error("Energy targeting exhausted eight corrections")
    current
end

function export_linear_reference(setup, target, prefix, provenance)
    p = setup.primary
    references, energies, residuals, factors = Dict(), Dict(), Dict(), Dict()
    tolerances, ulp_residuals = Dict(), Dict()
    for (name, modes, fundamental_m, level, step) in (
        ("linear_primary", filter(mode -> last(mode) == 1, p.modes), 1, 2, p.options.step),
        ("linear_secondary", setup.secondary_modes, 2, 3, 2p.options.step))
        (; inner, outer, match) = p.options
        # Locate the root, then measure its representable-frequency residual floor.
        # Unit outer amplitude makes the secondary's absolute residual sensitive
        # to a single ULP, as in the existing secondary continuation solver.
        coarse = solve_linear_mode(complex(p.alpha * (1 - p.alpha^2/(2level^2))),
            p.spin, p.alpha, first.(modes), fundamental_m, setup.couplings;
            step, inner, outer, match, tolerance=1e-7)
        residual_at(w) = matching_residual(w, coarse.horizon_amplitudes,
            coarse.outer_amplitudes, p.spin, p.alpha, first.(modes), fundamental_m,
            setup.couplings; step, inner, outer, match)
        center = residual_at(coarse.frequency)
        ulp = eps(real(coarse.frequency))
        floor = max(norm(residual_at(coarse.frequency+ulp)-center),
                    norm(residual_at(coarse.frequency-ulp)-center))
        tolerance = max(1e-8, 2floor)
        solution = solve_linear_mode(coarse.frequency, p.spin, p.alpha, first.(modes),
            fundamental_m, setup.couplings; horizon_guess=coarse.horizon_amplitudes,
            outer_tail_guess=coarse.outer_amplitudes[2:end],
            step, inner, outer, match, tolerance)
        @assert solution.residual_norm <= tolerance
        profiles = matching_residual(solution.frequency, solution.horizon_amplitudes,
            solution.outer_amplitudes, p.spin, p.alpha, first.(modes), fundamental_m,
            setup.couplings; step, inner, outer, match, save_profile=true)
        energy_of(radial) = sum(integrate_energy_branch(branch, solution.frequency,
            p.spin, p.alpha, modes, setup.polar, setup.azimuthal, setup.grid;
            fundamental_m, potential=:quadratic)
            for branch in (radial.from_horizon, radial.from_outer))
        energy = energy_of(profiles)
        @assert isfinite(energy) && energy > 0
        factor = sqrt(target / energy)
        scaled = (; (side => (; positions=getproperty(profiles, side).positions,
            states=factor .* getproperty(profiles, side).states)
            for side in (:from_horizon, :from_outer))...,
            residual=factor .* profiles.residual)
        measured = energy_of(scaled)
        @assert abs(measured / target - 1) < 1e-12
        references[name] = (; solution, profiles=scaled, modes)
        energies[name], residuals[name], factors[name] = measured, solution.residual_norm, factor
        tolerances[name], ulp_residuals[name] = tolerance, floor
        println("LINEAR ", name, " energy=", measured, " residual=", solution.residual_norm,
                " tolerance=", tolerance, " frequency_ulp_residual=", floor); flush(stdout)
    end
    metadata = Dict("passed"=>true, "normalization"=>"equal_cloud_energies",
        "alpha"=>p.alpha, "spin"=>p.spin, "inner"=>p.options.inner, "outer"=>p.options.outer,
        "reference_energy"=>target, "reference_curve_metadata"=>Dict("energies_by_curve"=>energies),
        "matching_residuals"=>residuals, "amplitude_factors"=>factors,
        "matching_tolerances"=>tolerances, "frequency_ulp_residuals"=>ulp_residuals,
        "spatial_convergence_established"=>false, "provenance"=>provenance)
    serialize(prefix * ".jls", (; references, metadata))
    open(prefix * ".csv", "w") do stream
        println(stream, "curve,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for name in ("linear_primary", "linear_secondary")
            reference = references[name]
            sample = fundamental_sampler(reference.profiles)
            ell, m = first(reference.modes)
            for radius in p.options.inner:0.25:p.options.outer
                state = sample(radius)
                @printf(stream, "%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    name, ell, m, radius, real(state.value), imag(state.value), abs(state.value),
                    real(state.slope), imag(state.slope))
            end
        end
    end
    open(prefix * ".toml", "w") do stream
        TOML.print(stream, metadata; sorted=true)
    end
end

function background_main()
    length(ARGS) in (2,3) || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY [--linear-only|--omiya-only]")
    length(ARGS) == 3 && ARGS[3] ∉ ("--linear-only", "--omiya-only") && error("Unknown option")
    linear_only = length(ARGS) == 3 && ARGS[3] == "--linear-only"
    omiya_only = length(ARGS) == 3 && ARGS[3] == "--omiya-only"
    inputs, output = abspath.(ARGS[1:2])
    ispath(output) && error("Preserve existing output; choose a fresh directory")
    config = TOML.parsefile(joinpath(inputs, "cases.toml"))
    cases = config["cases"]
    primary_path, secondary_path = joinpath(inputs, "primary.jls"), joinpath(inputs, "secondary.jls")
    saved_secondary = deserialize(secondary_path)
    seed_setup = load_primary_secondary_setup(primary_path; secondary_modes=saved_secondary.secondary_modes)
    seed = seed_setup.primary
    @assert Symbol(seed.potential) == :cosine && saved_secondary.reading == "averaged"
    @assert seed.spin == saved_secondary.spin && seed.alpha == saved_secondary.alpha
    @assert seed.options == saved_secondary.primary_options
    spin = only(unique(case["spin"] for case in cases))
    @assert isfinite(spin) && abs(spin) < 1
    for case in cases
        # This driver transports spin; other background changes need their own continuation.
        @assert case["alpha"] == seed.alpha
        @assert case["inner"] == seed.options.inner && case["outer"] == seed.options.outer
    end
    provenance = Dict("created_utc"=>string(now(UTC)), "slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
        "seed_spin"=>seed.spin,
        "input_sha256"=>Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in
            (primary_path, secondary_path, joinpath(inputs,"cases.toml"))))
    mkpath(output)
    targets = linear_only ? Float64[] : unique(case["target_energies"][1] for case in cases)
    for (index, target) in enumerate(targets)
        prefix = joinpath(output, "background$index")
        current = primary_target(seed_setup, target, prefix * "-seed")
        setup = seed_setup
        stages = max(1, ceil(Int, abs(spin-seed.spin)/0.02))
        for stage in 1:stages
            next_spin = stage == stages ? spin : seed.spin + (spin-seed.spin)*stage/stages
            setup = merge(seed_setup, (; primary=merge(seed, (; spin=next_spin))))
            if next_spin != seed.spin
                current = correct_reference_primary(setup, current, current.amplitude,
                    prefix * "-spin$stage")
                current = match_reference_energy(setup, current, target,
                    prefix * "-spin$stage-energy")
            end
        end
        serialize(prefix * "-primary.jls", (; current, spin, alpha=seed.alpha))
        secondary = secondary_target(setup, saved_secondary, current, prefix)
        for case in cases
            case["target_energies"][1] == target || continue
            export_reference(setup, current, secondary, case, joinpath(output, case["name"]), provenance)
        end
    end
    setup = merge(seed_setup, (; primary=merge(seed, (; spin))))
    if !omiya_only
        export_linear_reference(setup, config["reference_energy"], joinpath(output, "linear_reference"), provenance)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && background_main()
