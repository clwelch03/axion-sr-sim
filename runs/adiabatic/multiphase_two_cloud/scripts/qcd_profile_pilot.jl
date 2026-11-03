# One fixed-amplitude nonlinear solve, seeded by an energy-scaled linear mode.
# Usage: julia --project=. qcd_profile_pilot.jl CONFIG CASE_INDEX OUTPUT_DIRECTORY
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/single_cloud/load.jl"))

function pilot_metadata(path, data)
    open(path * ".tmp", "w") do stream
        TOML.print(stream, data; sorted=true)
    end
    mv(path * ".tmp", path; force=true)
end

function pilot_energy(profiles, frequency, modes, setup; potential)
    sum(integrate_energy_branch(branch, frequency, setup.spin, setup.alpha,
        modes, setup.polar, setup.azimuthal, setup.grid;
        fundamental_m=setup.fundamental_m, potential)
        for branch in (profiles.from_horizon, profiles.from_outer))
end

function pilot_flux(branch, frequency, setup)
    state = view(branch.states, :, 1)
    geometry = tortoise_geometry(first(branch.positions), setup.spin)
    count = length(setup.modes)
    4pi * geometry.radial_factor * sum(real(im * (mode[2] / setup.fundamental_m) *
        frequency * state[index] * conj(state[count+index]))
        for (index, mode) in enumerate(setup.modes))
end

function pilot_export_profiles(path, nonlinear, linear, modes, linear_modes)
    open(path, "w") do stream
        println(stream, "curve,branch,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for (curve, profiles, active_modes) in
                (("nonlinear", nonlinear, modes), ("linear", linear, linear_modes))
            for side in (:from_horizon, :from_outer)
                branch = getproperty(profiles, side)
                for (j, position) in enumerate(branch.positions), (i, (degree, m)) in enumerate(active_modes)
                    value, slope = branch.states[i,j], branch.states[length(active_modes)+i,j]
                    @printf(stream, "%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        curve, string(side), degree, m, position, real(value), imag(value),
                        abs(value), real(slope), imag(slope))
                end
            end
        end
    end
end

function qcd_pilot_main()
    length(ARGS) == 3 || error("usage: CONFIG CASE_INDEX OUTPUT_DIRECTORY")
    config_path, output = abspath(ARGS[1]), abspath(ARGS[3])
    config = TOML.parsefile(config_path)
    case = config["cases"][parse(Int, ARGS[2])]
    seed_path = joinpath(dirname(config_path), case["linear_seed"])
    seed = deserialize(seed_path)
    fundamental_m = case["dominant_m"]
    @assert seed.alpha == case["alpha"] && seed.spin == case["spin"]
    @assert seed.report["passed"] && first(seed.solution.outer_amplitudes) == 1
    @assert all(mode -> mode[2] == fundamental_m, seed.modes)
    @assert seed.report["normalization"] == "first outer coefficient = 1"
    modes = Tuple{Int,Int}[Tuple(mode) for mode in case["modes"]]
    @assert modes[1] == (fundamental_m, fundamental_m)
    @assert all(mode -> mode in modes, seed.modes)
    maximum_degree = maximum(first, modes)
    grid = angular_grid(case["polar_nodes"], case["azimuthal_nodes"])
    polar = build_polar_values(maximum_degree, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_degree, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_degree, polar, grid)
    setup = (; seed.alpha, seed.spin, fundamental_m, modes, grid, polar, azimuthal, couplings)
    options = seed.options
    target = case["target_code_energy"]
    @assert isfinite(target) && target > 0
    prefix = joinpath(output, case["name"])
    mkpath(output)
    any(startswith(name, case["name"] * suffix) for name in readdir(output)
        for suffix in (".", "-")) && error("Preserve existing output; use a fresh directory")
    provenance = Dict("created_utc"=>string(now(UTC)),
        "slurm_job_id"=>get(ENV, "SLURM_JOB_ID", "local"),
        "array_task_id"=>get(ENV, "SLURM_ARRAY_TASK_ID", "none"),
        "source_snapshot"=>normpath(joinpath(@__DIR__, "../../../..")),
        "config_sha256"=>bytes2hex(sha256(read(config_path))),
        "linear_seed_sha256"=>bytes2hex(sha256(read(seed_path))))
    pilot_metadata(prefix * "-request.toml", Dict("case"=>case, "provenance"=>provenance))

    # The same energy definition, evaluated with the pilot's angular tables,
    # must reproduce the archived free-field normalization for either m=1 or 2.
    linear_energy = pilot_energy(seed.profiles, seed.solution.frequency, seed.modes,
        setup; potential=:quadratic)
    @assert isapprox(linear_energy, seed.energy; rtol=1e-10)
    amplitude = sqrt(target / linear_energy)
    horizon = zeros(ComplexF64, length(modes))
    outer = zeros(ComplexF64, length(modes))
    for (index, mode) in enumerate(seed.modes)
        destination = only(findall(==(mode), modes))
        horizon[destination] = amplitude * seed.solution.horizon_amplitudes[index]
        outer[destination] = amplitude * seed.solution.outer_amplitudes[index]
    end
    tolerance = config["matching_tolerance"]
    serialize(prefix * "-input.jls", (; case, modes, options, amplitude,
        frequency_guess=seed.solution.frequency, horizon, outer, tolerance, provenance))
    println("DIRECT PILOT ", case["name"], " A=", amplitude,
        " target_energy=", target, " modes=", modes, " options=", options)
    flush(stdout)
    progress = (iteration, residual) -> begin
        println("NEWTON ", case["name"], " iteration=", iteration, " residual=", residual)
        flush(stdout)
    end
    timing = @timed solve_nonlinear_mode(seed.solution.frequency, amplitude,
        seed.spin, seed.alpha, modes, couplings, polar, azimuthal, grid;
        horizon_guess=horizon, outer_tail_guess=outer[2:end], fundamental_m,
        potential=:cosine, options..., tolerance, max_iterations=config["max_iterations"],
        jacobian_backend=:forwarddiff, progress_callback=progress)
    solution = timing.value
    serialize(prefix * "-solution.jls", (; solution, amplitude, case, modes, options, provenance))
    profiles = nonlinear_matching_residual(solution.frequency, solution.horizon_amplitudes,
        solution.outer_amplitudes, seed.spin, seed.alpha, modes, couplings, polar,
        azimuthal, grid; fundamental_m, potential=:cosine, options..., save_profile=true)
    matching_norm = norm(profiles.residual) / amplitude
    @assert matching_norm <= tolerance
    @assert all(branch -> all(isfinite, branch.states), (profiles.from_horizon, profiles.from_outer))
    energy = pilot_energy(profiles, solution.frequency, modes, setup; potential=:cosine)
    @assert isfinite(energy) && energy > 0
    energy_error = energy / target - 1
    within_energy_window = abs(energy_error) <= config["relative_energy_window"]

    # Only the linear comparison is rescaled. The nonlinear field is the solved
    # field, and its measured energy/occupation are reported even if off target.
    linear_factor = sqrt(energy / linear_energy)
    linear_profiles = (; (side => (; positions=getproperty(seed.profiles, side).positions,
        states=linear_factor .* getproperty(seed.profiles, side).states)
        for side in (:from_horizon, :from_outer))...,
        residual=linear_factor .* seed.profiles.residual)
    comparison_energy = pilot_energy(linear_profiles, seed.solution.frequency, seed.modes,
        setup; potential=:quadratic)
    @assert isapprox(comparison_energy, energy; rtol=1e-12)
    peak_field = 0.0
    for branch in (profiles.from_horizon, profiles.from_outer), index in 1:10:length(branch.positions)
        values = Dict(mode=>branch.states[i,index] for (i, mode) in enumerate(modes))
        peak_field = max(peak_field, maximum(abs, reconstruct_field(values, polar, azimuthal, grid)))
    end
    horizon_flux = -pilot_flux(profiles.from_horizon, solution.frequency, setup)
    outer_flux = pilot_flux(profiles.from_outer, solution.frequency, setup)
    measured_mass_fraction = energy * case["code_energy_to_mass_fraction"]
    measured_occupation = case["target_occupation"] * (energy / target) *
        real(seed.solution.frequency) / real(solution.frequency)
    report = Dict("passed"=>within_energy_window, "radial_solve_passed"=>true,
        "within_energy_window"=>within_energy_window, "potential"=>"cosine",
        "case"=>case, "target_energy"=>target, "measured_energy"=>energy,
        "relative_energy_error"=>energy_error, "relative_energy_window"=>config["relative_energy_window"],
        "measured_mass_fraction"=>measured_mass_fraction,
        "occupation_estimate_from_energy_over_frequency"=>measured_occupation,
        "amplitude"=>amplitude, "matching_norm"=>matching_norm, "matching_tolerance"=>tolerance,
        "frequency_real"=>real(solution.frequency), "frequency_imaginary"=>imag(solution.frequency),
        "linear_frequency_real"=>real(seed.solution.frequency),
        "linear_frequency_imaginary"=>imag(seed.solution.frequency),
        "linear_comparison_energy"=>comparison_energy, "linear_amplitude_factor"=>linear_factor,
        "sampled_peak_dimensionless_field"=>peak_field, "horizon_energy_flux"=>horizon_flux,
        "outer_energy_flux"=>outer_flux, "iterations"=>solution.iterations,
        "solve_seconds"=>timing.time, "inner"=>options.inner, "outer"=>options.outer,
        "match"=>options.match, "step"=>options.step, "horizon_order"=>3, "outer_order"=>7,
        "nonlinear_profile_rescaled"=>false, "time_evolution_performed"=>false,
        "spatial_convergence_established"=>false, "provenance"=>provenance)
    serialize(prefix * ".jls", (; solution, profiles, linear_profiles,
        linear_solution=seed.solution, linear_modes=seed.modes, modes, options, amplitude, report))
    pilot_export_profiles(prefix * ".csv", profiles, linear_profiles, modes, seed.modes)
    pilot_metadata(prefix * ".toml", report)
    println("FINISHED ", case["name"], " energy=", energy, " target=", target,
        " relative_error=", energy_error, " within_window=", within_energy_window)
    flush(stdout)
    within_energy_window || error("Profile saved, but outside requested energy window; inspect before retrying")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && qcd_pilot_main()
