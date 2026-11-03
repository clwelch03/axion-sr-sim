# Test whether moving the shooting match point resolves the failed 322 pilot.
# Retain the same physical target, cosine potential, modes, and radial spacing.
include(joinpath(@__DIR__, "qcd_profile_pilot.jl"))

function retry_main()
    length(ARGS) == 3 || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY MATCH_POINT")
    original_inputs, output = abspath(ARGS[1]), abspath(ARGS[2])
    match_point = parse(Float64,ARGS[3])
    ispath(output) && error("Preserve existing output; choose a fresh directory")
    inputs = joinpath(output,"inputs")
    mkpath(inputs)
    original_config = TOML.parsefile(joinpath(original_inputs,"pilots.toml"))
    case = original_config["cases"][2]
    @assert case["dominant_m"] == 2
    original_seed_path = joinpath(original_inputs,case["linear_seed"])
    seed = deserialize(original_seed_path)
    options = (; seed.options...,match=match_point)
    grid = angular_grid(16,31)
    polar = build_polar_values(6,grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(6,grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(6,polar,grid)
    solve_from(solution,tolerance) = solve_linear_mode(solution.frequency,seed.spin,seed.alpha,
        first.(seed.modes),2,couplings; options...,horizon_guess=solution.horizon_amplitudes,
        outer_tail_guess=solution.outer_amplitudes[2:end],tolerance,max_iterations=20)
    coarse = solve_from(seed.solution,1e-5)
    residual_at(frequency,solution;save_profile=false) = matching_residual(frequency,
        solution.horizon_amplitudes,solution.outer_amplitudes,seed.spin,seed.alpha,
        first.(seed.modes),2,couplings;options...,save_profile)
    center = residual_at(coarse.frequency,coarse)
    frequency_floor = maximum(norm(residual_at(coarse.frequency+sign*eps(real(coarse.frequency)),coarse)-center)
        for sign in (-1,1))
    tolerance = max(1e-8,2frequency_floor)
    solution = solve_from(coarse,tolerance)
    profiles = residual_at(solution.frequency,solution;save_profile=true)
    @assert norm(profiles.residual) <= tolerance
    energy = sum(integrate_energy_branch(branch,solution.frequency,seed.spin,seed.alpha,seed.modes,
        polar,azimuthal,grid;fundamental_m=2,potential=:quadratic)
        for branch in (profiles.from_horizon,profiles.from_outer))
    report = copy(seed.report)
    merge!(report,Dict("energy"=>energy,"match"=>match_point,
        "frequency_real"=>real(solution.frequency),"frequency_imaginary"=>imag(solution.frequency),
        "matching_norm"=>norm(profiles.residual),"acceptance_tolerance"=>tolerance,
        "frequency_ulp_residual"=>frequency_floor,"parent_seed_sha256"=>bytes2hex(sha256(read(original_seed_path)))))
    serialize(joinpath(inputs,case["linear_seed"]),merge(seed,(;solution,profiles,options,energy,report)))
    config = Dict("cases"=>[case],"matching_tolerance"=>original_config["matching_tolerance"],
        "max_iterations"=>original_config["max_iterations"],
        "relative_energy_window"=>original_config["relative_energy_window"])
    config_path = joinpath(inputs,"pilots.toml")
    open(config_path,"w") do stream
        TOML.print(stream,config;sorted=true)
    end
    println("RETRY changing only match point from ",seed.options.match," to ",match_point); flush(stdout)
    empty!(ARGS); append!(ARGS,[config_path,"1",output])
    qcd_pilot_main()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && retry_main()
