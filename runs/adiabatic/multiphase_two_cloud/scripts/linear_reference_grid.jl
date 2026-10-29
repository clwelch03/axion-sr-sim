# Two existing linear shooting solves per background; no nonlinear solve or
# QCD mass choice is needed for these dimensionless reference profiles.
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
const repository = normpath(joinpath(@__DIR__, "../../../.."))
include(joinpath(repository, "src/adiabatic/single_cloud/load.jl"))

function write_linear_metadata(path, data)
    open(path * ".tmp", "w") do stream
        TOML.print(stream, data; sorted=true)
    end
    mv(path * ".tmp", path; force=true)
end

# Same finite-boundary energy-current formula as energy_and_fluxes, evaluated
# on the already saved linear profiles to avoid a duplicate radial integration.
function linear_boundary_flux(branch, frequency, spin, mode_count)
    state = view(branch.states, :, 1)
    geometry = tortoise_geometry(first(branch.positions), spin)
    return 4pi * geometry.radial_factor * sum(
        real(im * frequency * state[index] * conj(state[mode_count+index]))
        for index in 1:mode_count)
end

function linear_reference_grid_main()
    length(ARGS) == 3 || error("usage: CONFIG_TOML CASE_INDEX OUTPUT_DIRECTORY")
    config_path, output_directory = abspath(ARGS[1]), abspath(ARGS[3])
    config = TOML.parsefile(config_path)
    case = config["cases"][parse(Int, ARGS[2])]
    alpha, spin = case["alpha"], case["spin"]
    isfinite(alpha) && alpha > 0 && isfinite(spin) && abs(spin) < 1 || error("Invalid background")
    prefix = joinpath(output_directory, case["name"])
    mkpath(output_directory)
    any(startswith(name, case["name"] * suffix) for name in readdir(output_directory)
        for suffix in (".", "-")) && error("Preserve existing output; choose a fresh directory")
    provenance = Dict("source_snapshot"=>repository,
        "config_sha256"=>bytes2hex(sha256(read(config_path))),
        "slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "array_task_id"=>get(ENV,"SLURM_ARRAY_TASK_ID","none"),
        "created_utc"=>string(now(UTC)))
    write_linear_metadata(prefix * "-request.toml", Dict("case"=>case, "provenance"=>provenance))
    grid = angular_grid(16, 31)
    polar = build_polar_values(6, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(6, polar, grid)
    cloud_reports = Dict{String,Any}()
    for (name, azimuthal_number, degrees, step) in (
            ("primary",1,[1,3,5],0.1), ("secondary",2,[2,4,6],0.2))
        cloud_prefix = prefix * "-" * name
        modes = [(degree,azimuthal_number) for degree in degrees]
        principal_number = azimuthal_number + 1
        frequency_guess = complex(alpha * (1 - alpha^2 / (2principal_number^2)))
        options = (; inner=-100.0, outer=Float64(case["outer"]), match=5.0, step)
        tolerance = 1e-8
        try
            println("SOLVING ",case["name"]," ",name," guess=",frequency_guess," options=",options)
            flush(stdout)
            timing = @timed solve_linear_mode(frequency_guess, spin, alpha,
                degrees, azimuthal_number, couplings; options..., tolerance, max_iterations=20)
            solution = timing.value
            serialize(cloud_prefix * "-solution.jls", (; solution, alpha, spin, modes, options))
            solution.residual_norm <= tolerance || error("Linear matching tolerance failed")
            profiles = matching_residual(solution.frequency, solution.horizon_amplitudes,
                solution.outer_amplitudes, spin, alpha, degrees, azimuthal_number, couplings;
                options..., save_profile=true)
            norm(profiles.residual) <= tolerance || error("Fresh radial reconstruction failed matching")
            all(branch -> all(isfinite, branch.states), (profiles.from_horizon,profiles.from_outer)) ||
                error("Nonfinite linear profile")
            energy = sum(integrate_energy_branch(branch, solution.frequency, spin, alpha,
                modes, polar, azimuthal, grid; fundamental_m=azimuthal_number, potential=:quadratic)
                for branch in (profiles.from_horizon,profiles.from_outer))
            isfinite(energy) && energy > 0 || error("Invalid free-field energy")
            horizon_flux = -linear_boundary_flux(profiles.from_horizon, solution.frequency, spin, length(modes))
            outer_flux = linear_boundary_flux(profiles.from_outer, solution.frequency, spin, length(modes))
            frequency_energy_rate = 2imag(solution.frequency) * energy
            balance_scale = max(abs(horizon_flux),abs(outer_flux),abs(frequency_energy_rate),eps(Float64))
            horizon_velocity = spin / (2kerr_horizons(spin).r_plus)
            report = Dict("passed"=>true,"alpha"=>alpha,"spin"=>spin,"m"=>azimuthal_number,
                "l_values"=>degrees,"principal_number_seed"=>principal_number,
                "frequency_real"=>real(solution.frequency),"frequency_imaginary"=>imag(solution.frequency),
                "frequency_guess_real"=>real(frequency_guess),"matching_norm"=>solution.residual_norm,
                "reconstructed_matching_norm"=>norm(profiles.residual),"requested_tolerance"=>tolerance,
                "iterations"=>solution.iterations,"solve_seconds"=>timing.time,
                "normalization"=>"first outer coefficient = 1",
                "energy"=>energy,"unit_energy_amplitude_factor"=>1/sqrt(energy),
                "horizon_energy_flux"=>horizon_flux,"outer_energy_flux"=>outer_flux,
                "relative_finite_domain_balance_error"=>(horizon_flux+outer_flux+frequency_energy_rate)/balance_scale,
                "superradiant_margin"=>azimuthal_number*horizon_velocity-real(solution.frequency),
                "inner"=>options.inner,"outer"=>options.outer,"match"=>options.match,"step"=>step,
                "horizon_order"=>3,"outer_order"=>7,"polar_nodes"=>16,"azimuthal_nodes"=>31,
                "spatial_convergence_established"=>false,"provenance"=>provenance)
            serialize(cloud_prefix * ".jls", (; solution, profiles, alpha, spin, modes, options,
                energy, horizon_flux, outer_flux, report))
            open(cloud_prefix * ".csv", "w") do stream
                println(stream,"branch,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
                for boundary in (:from_horizon,:from_outer)
                    branch = getproperty(profiles,boundary)
                    for (j,position) in enumerate(branch.positions), (i,degree) in enumerate(degrees)
                        field,slope = branch.states[i,j],branch.states[length(modes)+i,j]
                        @printf(stream,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                            string(boundary),degree,azimuthal_number,position,
                            real(field),imag(field),abs(field),real(slope),imag(slope))
                    end
                end
            end
            write_linear_metadata(cloud_prefix * ".toml", report)
            cloud_reports[name] = report
            println("ACCEPTED ",case["name"]," ",name," omega=",solution.frequency,
                " matching=",solution.residual_norm," E=",energy); flush(stdout)
        catch exception
            open(cloud_prefix * "-failure.txt","w") do stream
                showerror(stream,exception,catch_backtrace())
            end
            cloud_reports[name] = Dict("passed"=>false,"error"=>sprint(showerror,exception))
            println("FAILED ",case["name"]," ",name,": ",sprint(showerror,exception)); flush(stdout)
        end
    end
    passed = all(report["passed"] for report in values(cloud_reports))
    write_linear_metadata(prefix * ".toml", Dict("passed"=>passed,"alpha"=>alpha,
        "spin"=>spin,"clouds"=>cloud_reports,"provenance"=>provenance))
    passed || error("At least one linear mode failed; inspect saved per-cloud diagnostics")
end

if abspath(PROGRAM_FILE) == @__FILE__
    linear_reference_grid_main()
end
