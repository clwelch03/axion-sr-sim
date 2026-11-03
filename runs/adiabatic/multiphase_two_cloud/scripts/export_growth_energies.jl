# Reconstruct accepted profiles; do not evolve or rematch any trajectory.
# Reuse the per-cloud convention already used for energy-matched comparisons.
include("equal_energy_target.jl")

function growth_cloud_energies(radial, frequencies, problem, tables)
    primary_rows = vcat(tables.primary_indices, tables.mode_count .+ tables.primary_indices)
    secondary_rows = vcat(tables.secondary_indices, tables.mode_count .+ tables.secondary_indices)
    energies = zeros(2)
    for raw in (radial.from_horizon, radial.from_outer)
        branch = ordered_radial_branch(raw, tables.mode_count)
        densities = map(eachindex(branch.positions)) do index
            position = branch.positions[index]
            primary = branch.states[primary_rows, index]
            secondary = branch.states[secondary_rows, index]
            primary_values = Dict(mode => primary[i] for (i, mode) in enumerate(tables.primary_modes))
            e1 = energy_per_tortoise_radius(primary, position, frequencies[1],
                problem.spin, problem.alpha, tables.primary_modes, tables.polar,
                tables.azimuthal, tables.grid; potential=problem.potential)
            e2 = comparison_secondary_energy_density(secondary, position, primary_values,
                frequencies[2], problem.spin, problem.alpha, tables.secondary_modes,
                tables.polar, tables.azimuthal, tables.grid; potential=problem.potential)
            [e1, e2]
        end
        energies .+= [integrate_uniform(getindex.(densities, i), branch.spacing) for i in 1:2]
    end
    all(e -> isfinite(e) && e > 0, energies) || error("Invalid per-cloud energy")
    energies
end

function growth_linear_energies(problem, frequencies, amplitudes, tables)
    map(1:2) do cloud
        (; inner, outer, match, step, horizon_order, outer_order) = problem.integration
        linear = solve_linear_mode(frequencies[cloud], problem.spin, problem.alpha,
            [cloud], cloud, problem.couplings; inner, outer, match, step, tolerance=1e-9)
        isapprox(linear.frequency, frequencies[cloud]; rtol=1e-12) || error("Linear frequency changed")
        radial = nonlinear_matching_residual(linear.frequency,
            amplitudes[cloud] .* linear.horizon_amplitudes,
            ComplexF64[amplitudes[cloud]], problem.spin, problem.alpha,
            [(cloud, cloud)], problem.couplings, tables.polar, tables.azimuthal, tables.grid;
            inner, outer, match, step, horizon_order, outer_order,
            fundamental_m=cloud, potential=:quadratic, save_profile=true)
        sum(integrate_energy_branch(branch, real(linear.frequency), problem.spin,
            problem.alpha, [(cloud, cloud)], tables.polar, tables.azimuthal, tables.grid;
            fundamental_m=cloud, potential=:quadratic)
            for branch in (radial.from_horizon, radial.from_outer))
    end
end

function export_growth_energies(cosine, chiral, output; stride=1)
    stride >= 1 || error("Require a positive sample stride")
    any(ispath(output * suffix) for suffix in
        ("-cosine.csv", "-chiral.csv", "-postprocessing.toml")) && error("Choose a fresh output prefix")
    mkpath(dirname(output))
    input_hashes = Dict{String,String}()
    read_input(path) = begin
        bytes = read(path)
        input_hashes[abspath(path)] = bytes2hex(sha256(bytes))
        bytes
    end
    metadata = Dict{String,Any}("energy_unit" => "Fa^2 M", "time_unit" => "t/M",
        "convention" => "E211: isolated nonlinear primary; E322: quadratic secondary with phase-averaged U''(Phi211); mixed modes excluded",
        "carrier_frequency_convention" => "real frequencies; slow amplitude derivatives omitted, as in joint total energy",
        "normalization" => "shared outer-amplitude seeds, not equal initial energies",
        "linear_convention" => "free quadratic energies of linear profiles on the joint domain; E_i(t)=E_i(0)*exp(2*Im(omega_i)*t/M)",
        "stride" => stride, "created_utc" => string(now(UTC)))
    for (name, prefix) in (("cosine", cosine), ("chiral", chiral))
        report = TOML.parse(String(read_input(prefix * "-report.toml")))
        report["termination"] == "final_time" || error("Require a completed $name trajectory")
        problem, initial, frequencies = deserialize(IOBuffer(read_input(prefix * "-initial.jls")))
        steps = report["steps"]
        history = deserialize(IOBuffer(read_input(prefix * "-steps$steps.jls")))
        read_input(prefix * "-steps$steps.csv")
        problem.potential == Symbol(name) || error("Wrong potential")
        length(history.times) == steps + 1 && last(history.times) == report["requested_final_time"] || error("Incomplete time grid")
        tables = cloud_energy_tables(problem)
        # These growth runs have the fundamental primary plus two mixed modes.
        tables.primary_modes == [(1, 1)] && tables.secondary_modes == [(2, 2)] || error("Unexpected cloud support")
        cache = phase_averaged_observable_cache(problem)
        indices = sort(unique(vcat(collect(1:stride:length(history.times)), length(history.times))))
        check_indices = Set([first(indices), indices[cld(length(indices), 2)], last(indices)])
        checks = Dict{String,Any}[]
        output_csv = output * "-$name.csv"
        started = time()
        open(output_csv * ".tmp", "w") do stream
            println(stream, "step,time,A1,A2,E211,E322,total_energy,unassigned_energy")
            for index in indices
                sample = history.samples[index]
                profile = unpack_two_phase_profile(sample.solution.state, sample.solution.amplitudes, problem)
                radial = two_phase_matching_residual(profile, problem; save_radial_result=true)
                energies = growth_cloud_energies(radial, real.(profile.frequencies), problem, tables)
                total = sample.charges[1]
                println(stream, join([index-1; history.times[index]; sample.solution.amplitudes;
                    energies; total; total-sum(energies)], ','))
                if index in check_indices
                    reconstructed = phase_averaged_integrated_observables(radial,
                        real.(profile.frequencies), cache; potential=problem.potential).energy
                    total_error = abs(reconstructed/total - 1)
                    total_error < 1e-10 || error("Reconstructed total disagrees with checkpoint")
                    finer = comparison_cloud_energies(reference_sampler(radial),
                        real.(profile.frequencies), problem, tables; step=0.2)
                    quadrature_error = maximum(abs.(finer ./ energies .- 1))
                    quadrature_error < 1e-5 || error("Energy quadrature check failed")
                    push!(checks, Dict("step"=>index-1, "total_relative_error"=>total_error,
                        "quadrature_relative_difference"=>quadrature_error,
                        "cloud_energies"=>energies, "unassigned_energy"=>total-sum(energies)))
                end
                if index == first(indices) || (index-1) % 50 == 0 || index == last(indices)
                    println(name, " step=", index-1, "/", steps, " E=", energies,
                        " elapsed_s=", round(time()-started; digits=1)); flush(stdout); flush(stream)
                end
            end
        end
        mv(output_csv * ".tmp", output_csv)
        metadata[name] = Dict("prefix"=>abspath(prefix), "samples"=>length(indices),
            "potential"=>name, "alpha"=>problem.alpha, "spin"=>problem.spin,
            "spatial_settings"=>report["spatial_settings"], "checks"=>checks,
            "output_sha256"=>bytes2hex(sha256(read(output_csv))))
        if name == "cosine"
            linear_energies = growth_linear_energies(problem, frequencies, initial.solution.amplitudes, tables)
            maximum(abs.(linear_energies ./ first(checks)["cloud_energies"] .- 1)) < 1e-6 || error("Dilute energies disagree with free profiles")
            metadata["linear_initial_energies"] = linear_energies
            metadata["linear_amplitude_rates"] = imag.(frequencies)
        end
    end
    metadata["input_sha256"] = input_hashes
    source_paths = [@__FILE__, joinpath(@__DIR__, "equal_energy_target.jl"),
        joinpath(@__DIR__, "comparison_secondary_energy.jl"), joinpath(@__DIR__, "collocation_target.jl")]
    for library in ("core", "single_cloud", "primary_secondary", "multiphase")
        for (directory, _, files) in walkdir(joinpath(repository, "src/adiabatic", library))
            append!(source_paths, [joinpath(directory, file) for file in files if endswith(file, ".jl")])
        end
    end
    metadata["source_sha256"] = Dict(abspath(path)=>bytes2hex(sha256(read(path))) for path in source_paths)
    open(output * "-postprocessing.toml", "w") do stream
        TOML.print(stream, metadata; sorted=true)
    end
    println("Completed energy export: ", output); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (3, 4) || error("usage: COSINE_PREFIX CHIRAL_PREFIX OUTPUT_PREFIX [STRIDE]")
    export_growth_energies(abspath.(ARGS[1:3])...; stride=length(ARGS)==4 ? parse(Int, ARGS[4]) : 1)
end
