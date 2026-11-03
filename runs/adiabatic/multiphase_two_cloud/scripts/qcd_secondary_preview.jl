# Complete the two-cloud linear/Omiya preview on the accepted 211 background.
# Usage: julia --project=. qcd_secondary_preview.jl INPUT_DIRECTORY OUTPUT_PREFIX
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

function scaled_preview_profiles(profiles, factor)
    return (; (side => (; positions=getproperty(profiles, side).positions,
        states=factor .* getproperty(profiles, side).states)
        for side in (:from_horizon, :from_outer))...,
        residual=factor .* profiles.residual)
end

function preview_secondary_energy(profiles, frequency, primary, setup)
    energy = 0.0
    for side in (:from_horizon, :from_outer)
        branch, background = getproperty(profiles, side), getproperty(primary.profiles, side)
        order = sortperm(branch.positions)
        density = [secondary_energy_per_tortoise_radius(
            view(branch.states, :, i), branch.positions[i],
            primary_values_at(background, branch.positions[i], primary.modes),
            frequency, setup.spin, setup.alpha, setup.modes,
            setup.polar, setup.azimuthal, setup.grid) for i in order]
        energy += integrate_uniform(density, abs(branch.positions[2] - branch.positions[1]))
    end
    return energy
end

function main()
    length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_PREFIX")
    input_directory, prefix = abspath.(ARGS)
    any(ispath(prefix * suffix) for suffix in (".jls", ".toml", ".csv", "-secondary.jls")) &&
        error("Preserve existing results; choose a new prefix")
    mkpath(dirname(prefix))
    primary_path = joinpath(input_directory, "primary.jls")
    linear_path = joinpath(input_directory, "linear_secondary.jls")
    population_path = joinpath(input_directory, "population.toml")
    primary, linear = deserialize(primary_path), deserialize(linear_path)
    population = TOML.parsefile(population_path)
    @assert primary.report["passed"] && linear.report["passed"]
    case = primary.report["case"]
    dominant_mode = case["dominant_mode"]
    minor_mode = dominant_mode == "211" ? "322" : "211"
    minor_m = 3-case["dominant_m"]
    @assert case["alpha"] == linear.alpha && case["spin"] == linear.spin
    grid = angular_grid(case["polar_nodes"], case["azimuthal_nodes"])
    maximum_degree = max(maximum(first, primary.modes), maximum(first, linear.modes))
    polar = build_polar_values(maximum_degree, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_degree, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_degree, polar, grid)
    setup = (; linear.spin, linear.alpha, modes=linear.modes, polar, azimuthal, grid)
    solve_from(seed, tolerance) = solve_secondary_mode(seed.frequency,
        primary.profiles, primary.modes, linear.spin, linear.alpha, linear.modes,
        couplings, polar, azimuthal, grid; horizon_guess=seed.horizon_amplitudes,
        outer_tail_guess=seed.outer_amplitudes[2:end], tolerance,
        jacobian_backend=:forwarddiff, reading=:averaged)
    residual_at(frequency, solution; save_profile=false) = secondary_matching_residual(
        frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        primary.profiles, primary.modes, linear.spin, linear.alpha, linear.modes,
        couplings, polar, azimuthal, grid; save_profile, reading=:averaged)

    # Solve at unit secondary amplitude: the Omiya probe equation is linear in
    # this amplitude. Scaling afterward is exact for that equation and energy.
    println("Solving ",minor_mode," on the saved nonlinear ",dominant_mode," background"); flush(stdout)
    started = time()
    coarse = solve_from(linear.solution, 1e-5)
    center = residual_at(coarse.frequency, coarse)
    frequency_step = eps(real(coarse.frequency))
    frequency_floor = maximum(norm(residual_at(coarse.frequency + sign * frequency_step, coarse) - center)
        for sign in (-1, 1))
    tolerance = max(1e-8, 2frequency_floor)
    solution = solve_from(coarse, tolerance)
    unit_profiles = residual_at(solution.frequency, solution; save_profile=true)
    @assert norm(unit_profiles.residual) <= tolerance
    serialize(prefix * "-secondary.jls", (; solution, unit_profiles, tolerance, frequency_floor))
    unit_energy = preview_secondary_energy(unit_profiles, solution.frequency, primary, setup)
    @assert isfinite(unit_energy) && unit_energy > 0

    # Recover the finite microscopic reference recorded in the population
    # diagnostic. This is not a precise prediction of a classical occupation.
    target_energy = population["secondary_target_code_energy"]
    @assert isfinite(target_energy) && target_energy > 0
    omiya_factor = sqrt(target_energy / unit_energy)
    linear_factor = sqrt(target_energy / linear.energy)
    omiya_secondary = scaled_preview_profiles(unit_profiles, omiya_factor)
    linear_secondary = scaled_preview_profiles(linear.profiles, linear_factor)
    omiya_energy = preview_secondary_energy(omiya_secondary, solution.frequency, primary, setup)
    linear_energy = sum(integrate_energy_branch(branch, linear.solution.frequency,
        linear.spin, linear.alpha, linear.modes, polar, azimuthal, grid;
        fundamental_m=minor_m, potential=:quadratic)
        for branch in (linear_secondary.from_horizon, linear_secondary.from_outer))
    @assert isapprox(omiya_energy, target_energy; rtol=1e-12, atol=0)
    @assert isapprox(linear_energy, target_energy; rtol=1e-12, atol=0)
    @assert isapprox(primary.report["measured_energy"],
        primary.report["linear_comparison_energy"]; rtol=1e-12)
    references = (
        (; cloud=dominant_mode, model="linear", profiles=primary.linear_profiles, modes=primary.linear_modes),
        (; cloud=dominant_mode, model="omiya", profiles=primary.profiles, modes=primary.modes),
        (; cloud=minor_mode, model="linear", profiles=linear_secondary, modes=linear.modes),
        (; cloud=minor_mode, model="omiya", profiles=omiya_secondary, modes=linear.modes))
    open(prefix * ".csv", "w") do stream
        println(stream, "cloud,model,branch,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for reference in references, side in (:from_horizon, :from_outer)
            branch = getproperty(reference.profiles, side)
            @assert all(isfinite, branch.states)
            for (j, position) in enumerate(branch.positions), (i, (degree, m)) in enumerate(reference.modes)
                value, slope = branch.states[i,j], branch.states[length(reference.modes)+i,j]
                @printf(stream, "%s,%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    reference.cloud, reference.model, string(side), degree, m, position,
                    real(value), imag(value), abs(value), real(slope), imag(slope))
            end
        end
    end
    report = Dict("passed"=>true, "case"=>case, "population"=>population,
        "energy_"*dominant_mode=>primary.report["measured_energy"], "energy_"*minor_mode=>target_energy,
        "linear_energy_"*minor_mode=>linear_energy, "omiya_energy_"*minor_mode=>omiya_energy,
        "secondary_unit_energy"=>unit_energy, "secondary_unit_matching_norm"=>norm(unit_profiles.residual),
        "secondary_matching_tolerance"=>tolerance, "frequency_ulp_residual"=>frequency_floor,
        "secondary_frequency_real"=>real(solution.frequency),
        "secondary_frequency_imaginary"=>imag(solution.frequency),
        "secondary_omiya_amplitude_factor"=>omiya_factor,
        "secondary_linear_amplitude_factor"=>linear_factor,
        "joint_computed"=>false, "time_evolution_performed"=>false,
        "spatial_convergence_established"=>false, "solve_seconds"=>time()-started,
        "provenance"=>Dict("slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
            "created_utc"=>string(now(UTC)), "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
            "input_sha256"=>Dict(basename(path)=>bytes2hex(sha256(read(path)))
                for path in (primary_path, linear_path, population_path))))
    serialize(prefix * ".jls", (; references, secondary_solution=solution, report))
    open(prefix * ".toml", "w") do stream
        TOML.print(stream, report; sorted=true)
    end
    println("FINISHED both clouds, E211=", report["energy_211"], " E322=", report["energy_322"],
        " secondary matching=", norm(unit_profiles.residual)); flush(stdout)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
