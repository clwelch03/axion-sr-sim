# Validate actual campaign inputs and independent energy targeting on a small
# four-component cosine problem before starting the expensive parameter grid.
using Test
include("parameter_energy_target.jl")
length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY")
inputs, output = abspath.(ARGS)
mkpath(output)
seed = deserialize(joinpath(inputs, "initial.jls"))
config = TOML.parsefile(joinpath(inputs, "sweep.toml"))

@testset "Requested grid and independent coordinates" begin
    @test Set(case["alpha"] for case in config["cases"]) == Set([0.2, 0.3, 0.4])
    @test Set(case["spin"] for case in config["cases"]) == Set([0.8, 0.9, 0.95, 0.99])
    @test length(config["cases"]) == 60
    @test seed.assessment.passed
    tables = cloud_energy_tables(seed.problem)
    # The production primary includes m=5, so the corrected helper must return
    # exactly the same angular table as the running frozen campaign.
    @test tables.azimuthal == build_azimuthal_values(maximum(last, tables.primary_modes), tables.grid.azimuth_nodes)
    energies = measured_cloud_energies(seed; step=0.5)
    @test maximum(abs.(energies ./ seed.energies .- 1)) < 1e-5
    @test maximum(abs.(energies ./ config["reference_energy"] .- 1)) < 3e-5
    expected_pairs = Set([(1.0,1.0), (0.5,1.0), (1.0,0.5), (1.5,1.0), (1.0,1.5)])
    for alpha in (0.2,0.3,0.4), spin in (0.8,0.9,0.95,0.99)
        cases = filter(case -> case["alpha"] == alpha && case["spin"] == spin, config["cases"])
        @test Set(Tuple(case["energy_factors"]) for case in cases) == expected_pairs
        for case in cases
            validate_parameter_case(case)
            @test case["energies"] == config["reference_energy"] .* case["energy_factors"]
        end
        case = first(cases)
        problem = changed_background(seed.problem; alpha, spin, outer=case["outer"])
        prediction = energy_predictor(seed, seed.amplitudes; problem)
        @test problem.alpha == alpha && problem.spin == spin
        @test problem.source_basis === seed.problem.source_basis
        @test two_phase_profile_admissible(pack_two_phase_profile(prediction.profile, problem), problem)
        @test all(isfinite, two_phase_boundary_state(prediction.profile, problem, :horizon))
        @test all(isfinite, two_phase_boundary_state(prediction.profile, problem, :outer))
        @test all(position -> all(isfinite, prediction.guess(position)),
            range(problem.integration.inner, problem.integration.outer; length=31))
    end
    @test seed.problem.alpha == 0.3 && seed.problem.spin == 0.99
    @test seed.problem.integration.outer == 300.0
    invalid = merge(first(config["cases"]), Dict("energies"=>[1.0,0.0]))
    @test_throws ErrorException validate_parameter_case(invalid)
end

@testset "Unequal energy targets with a changed background" begin
    grid = angular_grid(8, 9)
    polar_values = build_polar_values(3, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(3, polar_values, grid)
    cloud_m, spin, alpha = [1,2], 0.99, 0.42
    components = [fourier_component(q, [ell], cloud_m) for (ell,q) in
        [(1,(1,0)), (2,(0,1)), (0,(2,-1)), (3,(-1,2))]]
    integration = (; inner=-30.0, outer=80.0, match=5.0, step=0.4, horizon_order=3, outer_order=7)
    options = (; inner=integration.inner, outer=integration.outer, match=integration.match, step=integration.step)
    first_mode = solve_linear_mode(0.4088+1e-7im, spin, alpha, [1], 1, couplings; options..., tolerance=1e-9)
    second_mode = solve_linear_mode(0.4159+1e-9im, spin, alpha, [2], 2, couplings; options..., tolerance=1e-9)
    amplitudes = [1e-4,1e-5]
    profile = TwoPhaseProfile([first_mode.frequency,second_mode.frequency],
        ComplexF64[amplitudes[1]*first_mode.horizon_amplitudes[1], amplitudes[2]*second_mode.horizon_amplitudes[1],0,0],
        ComplexF64[amplitudes...,0,0])
    problem = TwoPhaseProblem(; spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes=grid.azimuth_nodes, phase2_nodes=grid.azimuth_nodes, potential=:cosine,
        coefficient_scale=minimum(amplitudes), free_outer_coefficient_indices=[3,4],
        initial_bound=BitVector([true,true,true,false]), integration)
    radial = two_phase_matching_residual(profile, problem; save_radial_result=true)
    guess = reference_sampler(radial)
    scales = two_phase_collocation_scales(problem, amplitudes, profile, guess)
    result = solve_two_phase_collocation(problem, amplitudes, profile, guess;
        scales..., initial_step=1.0, tolerance=1e-5, nonlinear_tolerance=1e-7,
        ad_chunksize=8, time_limit=900.0)
    assessment = assess_two_phase_collocation(result, problem)
    @test assessment.passed
    current = energy_record(result, assessment, problem, 0.0)
    @test haskey(cloud_energy_tables(problem).azimuthal, 2)
    original_energies = measured_cloud_energies(current)
    targets = original_energies .* [0.7,1.3]
    changed = changed_background(problem; alpha=0.418, spin=0.98, outer=82.0)
    history = Dict{String,Any}[]
    prefix = joinpath(output,"unequal_energy_check")
    current = solve_parameter_step(current, changed, amplitudes, prefix, history)
    current = match_energy_pair(current, targets, prefix, history; final=true)
    current = solve_parameter_step(current, changed, current.amplitudes, prefix, history; tolerance=1e-6)
    current = match_energy_pair(current, targets, prefix, history;
        final=true, tolerance=1e-6, quadrature_step=0.5)
    measured = measured_cloud_energies(current; step=0.5)
    @test maximum(abs.(measured ./ targets .- 1)) < 1e-4
    @test measured[1] < original_energies[1] && measured[2] > original_energies[2]
    @test current.problem.alpha == 0.418 && current.problem.spin == 0.98
    @test current.problem.integration.outer == 82.0
    @test current.assessment.passed
    @test current.settings.tolerance == 1e-6
    @test maximum(abs.(current.assessment.states[3:4,:])) > 0
    export_parameter_case(current, targets, prefix, history, Dict("purpose"=>"small cosine validation"))
end
save_toml(joinpath(output,"validation.toml"), Dict("passed"=>true,
    "production_cases_checked"=>length(config["cases"]), "created_utc"=>string(now(UTC))))
println("PARAMETER/ENERGY VALIDATION PASSED"); flush(stdout)
