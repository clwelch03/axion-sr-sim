# Run the complete linear/Omiya/joint-probe sequence for one frozen QCD case.
include(joinpath(@__DIR__, "qcd_profile_pilot.jl"))
include(joinpath(@__DIR__, "qcd_joint_probe.jl"))

# Rematch a free reference when the dominant-322 setup changes its match point
# or the probe's RK4 spacing. The microscopic population is unchanged.
function rematch_comparison_seed(seed,options)
    grid = angular_grid(16,31)
    polar = build_polar_values(6,grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6,grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(6,polar,grid)
    fundamental_m = last(first(seed.modes))
    solve_from(solution,tolerance) = solve_linear_mode(solution.frequency,seed.spin,seed.alpha,
        first.(seed.modes),fundamental_m,couplings; options...,horizon_guess=solution.horizon_amplitudes,
        outer_tail_guess=solution.outer_amplitudes[2:end],tolerance,max_iterations=20)
    coarse = solve_from(seed.solution,1e-5)
    residual_at(frequency,solution;save_profile=false) = matching_residual(frequency,
        solution.horizon_amplitudes,solution.outer_amplitudes,seed.spin,seed.alpha,
        first.(seed.modes),fundamental_m,couplings;options...,save_profile)
    center = residual_at(coarse.frequency,coarse)
    frequency_floor = maximum(norm(residual_at(coarse.frequency+sign*eps(real(coarse.frequency)),coarse)-center)
        for sign in (-1,1))
    tolerance = max(1e-8,2frequency_floor)
    solution = solve_from(coarse,tolerance)
    profiles = residual_at(solution.frequency,solution;save_profile=true)
    @assert norm(profiles.residual) <= tolerance
    energy = sum(integrate_energy_branch(branch,solution.frequency,seed.spin,seed.alpha,seed.modes,
        polar,azimuthal,grid;fundamental_m,potential=:quadratic)
        for branch in (profiles.from_horizon,profiles.from_outer))
    report = copy(seed.report)
    merge!(report,Dict("energy"=>energy,"match"=>options.match,"step"=>options.step,
        "frequency_real"=>real(solution.frequency),"frequency_imaginary"=>imag(solution.frequency),
        "matching_norm"=>norm(profiles.residual),"acceptance_tolerance"=>tolerance,
        "frequency_ulp_residual"=>frequency_floor))
    return merge(seed,(;solution,profiles,options,energy,report))
end

function qcd_211_comparison_main()
    length(ARGS) == 3 || error("usage: CONFIG CASE_INDEX OUTPUT_DIRECTORY")
    config_path, index, output = abspath(ARGS[1]), parse(Int,ARGS[2]), abspath(ARGS[3])
    config = TOML.parsefile(config_path)
    case = config["cases"][index]
    case_output = joinpath(output,case["name"])
    ispath(case_output) && error("Preserve existing case output")
    mkpath(case_output)
    inputs = joinpath(case_output,"inputs")
    mkpath(inputs)
    pilot_config_path,pilot_index = config_path,index
    if case["dominant_m"] == 2
        seed = deserialize(joinpath(dirname(config_path),case["linear_seed"]))
        # 24 is close to the successful match=25 test, and makes both branch
        # lengths multiples of .4 so every probe half-step is a primary sample.
        rematched = rematch_comparison_seed(seed,(;seed.options...,match=24.0))
        serialize(joinpath(inputs,case["linear_seed"]),rematched)
        pilot_config = copy(config)
        pilot_config["cases"] = [case]
        pilot_config_path,pilot_index = joinpath(inputs,"pilot-case.toml"),1
        open(pilot_config_path,"w") do stream
            TOML.print(stream,pilot_config;sorted=true)
        end
    end
    empty!(ARGS); append!(ARGS,[pilot_config_path,string(pilot_index),case_output])
    qcd_pilot_main()
    primary_path = joinpath(case_output,case["name"]*".jls")
    primary = deserialize(primary_path)
    cp(primary_path,joinpath(inputs,"primary.jls"))
    secondary_path = joinpath(dirname(config_path),case["secondary_seed"])
    linear = deserialize(secondary_path)
    if case["dominant_m"] == 2
        linear = rematch_comparison_seed(linear,(;primary.options...,step=2primary.options.step))
    end
    if !hasproperty(linear,:solution)
        # The recovered alpha=.3/spin=.95 reference predates the grid archive
        # format. Reconstruct its unit-amplitude samples from its accepted root.
        solution = linear.references["linear_secondary"].solution
        @assert first(solution.outer_amplitudes) == 1
        grid = angular_grid(16,31)
        polar = build_polar_values(6,grid.cosine_theta_nodes)
        azimuthal = build_azimuthal_values(6,grid.azimuth_nodes)
        couplings = build_cosine_squared_couplings(6,polar,grid)
        modes = [(2,2),(4,2),(6,2)]
        options = (; primary.options...,step=0.2)
        profiles = matching_residual(solution.frequency,solution.horizon_amplitudes,
            solution.outer_amplitudes,case["spin"],case["alpha"],[2,4,6],2,couplings;
            options...,save_profile=true)
        @assert norm(profiles.residual) <= 1.0964964352532414e-7
        energy = sum(integrate_energy_branch(branch,solution.frequency,case["spin"],case["alpha"],
            modes,polar,azimuthal,grid;fundamental_m=2,potential=:quadratic)
            for branch in (profiles.from_horizon,profiles.from_outer))
        linear = (; solution,profiles,modes,options,energy,alpha=case["alpha"],spin=case["spin"],
            report=Dict("passed"=>true))
    end
    serialize(joinpath(inputs,"linear_secondary.jls"),linear)
    occupation = 10.0^case["log10_minor_occupation"]
    target_energy = occupation*real(linear.solution.frequency)/case["occupation_scale"] /
        case["code_energy_to_mass_fraction"]
    population = Dict("secondary_reference_occupation"=>occupation,
        "secondary_target_code_energy"=>target_energy,"seed_occupation"=>1.0,
        "original_log10_minor_occupation"=>case["log10_minor_occupation"],
        "occupation_scale"=>case["occupation_scale"],
        "code_energy_to_mass_fraction"=>case["code_energy_to_mass_fraction"],
        "conversion_frequency"=>real(linear.solution.frequency),
        "population_diagnostic_quartic_coefficient"=>0.3,
        "profile_potential"=>"cosine","profile_quartic_coefficient"=>1.0,
        "population_interpretation"=>"Preserved original first-peak seed diagnostic; microscopic occupation is illustrative, not a precise classical prediction.")
    open(joinpath(inputs,"population.toml"),"w") do stream
        TOML.print(stream,population;sorted=true)
    end
    empty!(ARGS); append!(ARGS,[inputs,joinpath(case_output,"preview")])
    main() # qcd_secondary_preview.jl
    cp(joinpath(case_output,"preview.jls"),joinpath(inputs,"preview.jls"))
    empty!(ARGS); append!(ARGS,[inputs,joinpath(case_output,"comparison")])
    joint_probe_main()
    println("ALL THREE MODELS FINISHED ",case["name"]); flush(stdout)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && qcd_211_comparison_main()
