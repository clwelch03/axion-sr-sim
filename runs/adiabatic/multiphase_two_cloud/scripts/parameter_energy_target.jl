# Continue the joint collocation solution in (alpha, spin, E1, E2).
# The two comparison energies are independent coordinates; the full interacting
# energy and angular momentum are measured separately at the accepted endpoint.
include("equal_energy_target.jl")

function changed_background(problem; alpha=problem.alpha, spin=problem.spin,
        outer=problem.integration.outer)
    fields = (; (name => getfield(problem, name) for name in fieldnames(typeof(problem)))...)
    integration = merge(problem.integration, (; outer=Float64(outer)))
    return TwoPhaseProblem(; merge(fields, (; alpha=Float64(alpha), spin=Float64(spin), integration))...)
end

function continuation_sampler(current)
    assessment, problem = current.assessment, current.problem
    return position -> begin
        if position <= last(assessment.positions)
            hermite_state(assessment.positions, assessment.states, position)
        else
            # Extend only the initial guess with the previous solution's outer
            # series. Every extended domain is subsequently solved and assessed.
            extended = changed_background(problem; outer=position)
            two_phase_boundary_state(assessment.profile, extended, :outer)
        end
    end
end

function energy_predictor(current, amplitudes; problem=current.problem)
    ratios = amplitudes ./ current.amplitudes
    factors = [prod(ratios .^ abs.(collect(q))) for (_, q) in problem.source_basis.mode_labels]
    sampler = continuation_sampler(current)
    guess = position -> sampler(position) .* vcat(factors, factors)
    old_profile = current.assessment.profile
    # Preserve omega/alpha only as a predictor, then solve both complex frequencies.
    frequencies = old_profile.frequencies .* (problem.alpha / current.problem.alpha)
    profile = TwoPhaseProfile(frequencies, old_profile.horizon_coefficients .* factors,
        old_profile.outer_coefficients .* factors)
    return (; guess, profile)
end

function save_toml(path, contents)
    open(path * ".tmp", "w") do stream
        TOML.print(stream, contents; sorted=true)
    end
    mv(path * ".tmp", path; force=true)
end

function measured_cloud_energies(current; step=1.0)
    return comparison_cloud_energies(continuation_sampler(current),
        current.assessment.profile.frequencies, current.problem,
        cloud_energy_tables(current.problem); step)
end

function solve_parameter_step(current, problem, amplitudes, prefix, history;
        tolerance=1e-5, maximum_solves=64)
    solve_number = length(history) + 1
    solve_number <= maximum_solves || error("Reached the bounded solve count: $maximum_solves")
    (; guess, profile) = energy_predictor(current, amplitudes; problem)
    two_phase_profile_admissible(pack_two_phase_profile(profile, problem), problem) ||
        error("Parameter predictor crossed a frequency branch")
    scales = two_phase_collocation_scales(problem, amplitudes, profile, guess)
    initial_step = tolerance <= 1e-6 ? 0.5 : 1.0
    maximum_intervals = max(1600, ceil(Int, 3 * (problem.integration.outer - problem.integration.inner)))
    println("SOLVE ", solve_number, " alpha=", problem.alpha, " spin=", problem.spin,
        " amplitudes=", amplitudes, " outer=", problem.integration.outer, " tolerance=", tolerance)
    flush(stdout)
    timing = @timed solve_two_phase_collocation(problem, amplitudes, profile, guess;
        scales..., initial_step, tolerance, nonlinear_tolerance=tolerance/100,
        max_mesh_intervals=maximum_intervals, ad_chunksize=min(28, length(scales.state_scales)),
        time_limit=1800.0, trace=true)
    result = timing.value
    assessment = assess_two_phase_collocation(result, problem)
    record = energy_record(result, assessment, problem, timing.time)
    serialize(prefix * "-solve$(solve_number).jls", record)
    row = Dict("solve"=>solve_number, "alpha"=>problem.alpha, "spin"=>problem.spin,
        "amplitudes"=>amplitudes, "outer"=>problem.integration.outer,
        "tolerance"=>tolerance, "passed"=>assessment.passed,
        "elapsed_seconds"=>timing.time, "mesh_points"=>length(record.mesh),
        "boundary_residual"=>assessment.boundary_residual,
        "defect"=>maximum(assessment.defect_by_state))
    push!(history, row)
    save_toml(prefix * "-history.toml", Dict("solves"=>history))
    assessment.passed || error("Radial assessment failed at solve $solve_number; checkpoint preserved")
    energies = measured_cloud_energies(record)
    all(value -> isfinite(value) && value > 0, energies) || error("Nonpositive cloud energy")
    row["energies"] = energies
    save_toml(prefix * "-history.toml", Dict("solves"=>history))
    return record
end

function match_energy_pair(current, targets, prefix, history; final=false,
        tolerance=1e-5, quadrature_step=1.0)
    jacobian = 2.0 .* Matrix{Float64}(I, 2, 2)
    previous_x = previous_error = nothing
    threshold = final ? 3e-5 : 3e-4
    for iteration in 0:8
        energies = measured_cloud_energies(current; step=quadrature_step)
        all(value -> isfinite(value) && value > 0, energies) || error("Nonpositive cloud energy")
        error_vector = log.(energies ./ targets)
        println("ENERGY targets=", targets, " measured=", energies,
            " relative_error=", energies ./ targets .- 1); flush(stdout)
        maximum(abs, error_vector) < threshold && return current
        iteration < 8 || error("Energy matching exceeded eight corrections")
        x = log.(current.amplitudes)
        if previous_x !== nothing
            dx, dy = x - previous_x, error_vector - previous_error
            dot(dx, dx) > 1e-20 && (jacobian += (dy - jacobian * dx) * transpose(dx) / dot(dx, dx))
        end
        correction = -(jacobian \ error_vector)
        all(isfinite, correction) || error("Nonfinite energy correction")
        correction *= min(1.0, 0.2 / max(maximum(abs, correction), eps()))
        previous_x, previous_error = x, error_vector
        current = solve_parameter_step(current, current.problem, exp.(x + correction), prefix, history; tolerance)
    end
    error("Energy matching did not return an accepted solution")
end

function validate_parameter_case(case)
    all(key -> haskey(case, key), ("name", "alpha", "spin", "energies", "outer")) ||
        error("Each case requires name, alpha, spin, energies, and outer")
    isfinite(case["alpha"]) && case["alpha"] > 0 || error("Require positive finite alpha")
    isfinite(case["spin"]) && abs(case["spin"]) < 1 || error("Require |spin| < 1")
    length(case["energies"]) == 2 && all(value -> isfinite(value) && value > 0, case["energies"]) ||
        error("Require two independent positive finite cloud energies")
    isfinite(case["outer"]) && case["outer"] >= 300 && isinteger(case["outer"]) ||
        error("Require an integer outer endpoint of at least 300")
    basename(case["name"]) == case["name"] && !(case["name"] in ("", ".", "..")) ||
        error("Require a simple case name")
end

function export_parameter_case(current, targets, prefix, history, provenance)
    problem, assessment = current.problem, current.assessment
    energies = measured_cloud_energies(current)
    fine_energies = measured_cloud_energies(current; step=0.5)
    maximum(abs.(fine_energies ./ energies .- 1)) < 5e-6 || error("Energy quadrature check failed")
    maximum(abs.(fine_energies ./ targets .- 1)) < 1e-4 || error("Refined energies missed their targets")
    indices = 1:4:length(assessment.positions)
    branch = (; positions=assessment.positions[indices], states=assessment.states[:, indices], spacing=1.0)
    total = integrate_two_phase_branch_observables(branch, problem.source_basis.mode_labels,
        problem.spin, problem.alpha, real.(assessment.profile.frequencies), problem.cloud_m,
        problem.polar_values, problem.grid, problem.phase1_nodes, problem.phase2_nodes;
        potential=problem.potential)
    record = merge(current, (; energies=fine_energies, target_energies=targets,
        history, total_observables=total, provenance))
    serialize(prefix * ".jls", record)
    metadata = Dict("passed"=>true, "alpha"=>problem.alpha, "spin"=>problem.spin,
        "target_energies"=>targets, "joint_cloud_energies"=>fine_energies,
        "relative_energy_errors"=>fine_energies ./ targets .- 1,
        "energy_convention"=>"isolated cosine primary; quadratic secondary on averaged primary; mixed modes excluded from labels",
        "joint_total_energy"=>total.energy, "joint_total_angular_momentum"=>total.angular_momentum,
        "horizon_energy_flux"=>-total.lower.radial_energy_flux,
        "outer_energy_flux"=>total.upper.radial_energy_flux,
        "horizon_angular_momentum_flux"=>-total.lower.radial_angular_momentum_flux,
        "outer_angular_momentum_flux"=>total.upper.radial_angular_momentum_flux,
        "cloud_amplitudes"=>current.amplitudes,
        "frequencies_real"=>real.(assessment.profile.frequencies),
        "frequencies_imaginary"=>imag.(assessment.profile.frequencies),
        "inner"=>problem.integration.inner, "outer"=>problem.integration.outer,
        "horizon_order"=>problem.integration.horizon_order, "outer_order"=>problem.integration.outer_order,
        "mode_labels"=>string.(problem.source_basis.mode_labels),
        "polar_nodes"=>length(problem.grid.cosine_theta_nodes),
        "phase_nodes"=>[length(problem.phase1_nodes), length(problem.phase2_nodes)],
        "mesh_points"=>length(current.mesh), "requested_tolerance"=>current.settings.tolerance,
        "boundary_residual"=>assessment.boundary_residual,
        "independent_off_mesh_defect"=>maximum(assessment.defect_by_state),
        "parameter_drift"=>assessment.parameter_drift,
        "spatial_convergence_established"=>false, "provenance"=>provenance)
    save_toml(prefix * ".toml", metadata)
    open(prefix * ".csv", "w") do stream
        println(stream, "ell,q1,q2,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        labels = problem.source_basis.mode_labels
        for (j, position) in enumerate(assessment.positions), (i, (ell, q)) in enumerate(labels)
            field, slope = assessment.states[i,j], assessment.states[length(labels)+i,j]
            @printf(stream, "%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                ell,q[1],q[2],position,real(field),imag(field),abs(field),real(slope),imag(slope))
        end
    end
    println("ACCEPTED alpha=", problem.alpha, " spin=", problem.spin,
        " energies=", fine_energies, " targets=", targets); flush(stdout)
end

function parameter_energy_main()
    length(ARGS) in (4, 5) || error("usage: SEED_JLS CONFIG_TOML CASE_INDEX OUTPUT_DIRECTORY [--check]")
    checking = length(ARGS) == 5
    checking && ARGS[5] != "--check" && error("Unknown option")
    seed_path, config_path = abspath.(ARGS[1:2])
    config = TOML.parsefile(config_path)
    case = config["cases"][parse(Int, ARGS[3])]
    validate_parameter_case(case)
    seed = deserialize(seed_path)
    seed.assessment.passed || error("Seed did not pass its radial assessment")
    seed_energies = measured_cloud_energies(seed; step=0.5)
    maximum(abs.(seed_energies ./ seed.energies .- 1)) < 1e-5 || error("Saved seed energy mismatch")
    targets = Float64.(case["energies"])
    steps = max(1, ceil(Int, abs(case["alpha"] - seed.problem.alpha)/0.02),
        ceil(Int, abs(case["spin"] - seed.problem.spin)/0.025),
        ceil(Int, maximum(abs, log.(targets ./ seed_energies))/0.25))
    println("CASE ", case["name"], " steps=", steps, " targets=", targets,
        " seed_energies=", seed_energies); flush(stdout)
    if checking
        target_problem = changed_background(seed.problem; alpha=case["alpha"],
            spin=case["spin"], outer=case["outer"])
        prediction = energy_predictor(seed, seed.amplitudes; problem=target_problem)
        @assert two_phase_profile_admissible(pack_two_phase_profile(prediction.profile, target_problem), target_problem)
        for side in (:horizon, :outer)
            @assert all(isfinite, two_phase_boundary_state(prediction.profile, target_problem, side))
        end
        for position in range(target_problem.integration.inner, target_problem.integration.outer; length=31)
            @assert all(isfinite, prediction.guess(position))
        end
        println("SETUP CHECK PASSED; no nonlinear solve performed")
        return
    end
    mkpath(abspath(ARGS[4]))
    prefix = joinpath(abspath(ARGS[4]), case["name"])
    any(startswith(name, case["name"] * suffix) for name in readdir(dirname(prefix))
        for suffix in (".", "-")) && error("Preserve existing case output; choose a fresh directory")
    provenance = Dict("seed_path"=>seed_path, "seed_sha256"=>bytes2hex(sha256(read(seed_path))),
        "config_sha256"=>bytes2hex(sha256(read(config_path))), "source_snapshot"=>repository,
        "slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "slurm_array_task_id"=>get(ENV,"SLURM_ARRAY_TASK_ID","none"), "created_utc"=>string(now(UTC)))
    save_toml(prefix * "-request.toml", Dict("case"=>case, "provenance"=>provenance, "stages"=>steps))
    history = Dict{String,Any}[]
    current = seed
    try
        if case["alpha"] == seed.problem.alpha && case["spin"] == seed.problem.spin &&
                case["outer"] == seed.problem.integration.outer &&
                maximum(abs.(seed_energies ./ targets .- 1)) < 3e-5 && seed.settings.tolerance <= 1e-6
            export_parameter_case(seed, targets, prefix, history, provenance)
            return
        end
        for stage in 1:steps
            fraction = stage/steps
            alpha = stage == steps ? case["alpha"] : seed.problem.alpha + fraction*(case["alpha"]-seed.problem.alpha)
            spin = stage == steps ? case["spin"] : seed.problem.spin + fraction*(case["spin"]-seed.problem.spin)
            # Increase the domain gradually; all endpoints remain compatible
            # with the 1, 0.5, and 0.25 quadrature/sample spacings.
            outer = ceil(seed.problem.integration.outer + fraction*(case["outer"]-seed.problem.integration.outer))
            problem = changed_background(current.problem; alpha, spin, outer)
            stage_targets = stage == steps ? targets : seed_energies .* exp.(fraction .* log.(targets ./ seed_energies))
            # changed_background only changes these settings. Reuse an accepted
            # solution when they are identical; match_energy_pair performs the
            # next actual amplitude change and its required radial correction.
            if problem.alpha != current.problem.alpha || problem.spin != current.problem.spin ||
                    problem.integration != current.problem.integration
                current = solve_parameter_step(current, problem, current.amplitudes, prefix, history)
            else
                println("Reusing accepted radial solution at unchanged background for stage ", stage)
                flush(stdout)
            end
            current = match_energy_pair(current, stage_targets, prefix, history; final=stage==steps)
            serialize(prefix * "-stage$(stage).jls", merge(current, (; energies=measured_cloud_energies(current))))
        end
        current = solve_parameter_step(current, current.problem, current.amplitudes, prefix, history; tolerance=1e-6)
        # A successful radial refinement may shift an observable more than its
        # targeting tolerance. Correct the energies using the final mesh tolerance
        # and quadrature, rather than accepting the coarser amplitudes unchanged.
        current = match_energy_pair(current, targets, prefix, history;
            final=true, tolerance=1e-6, quadrature_step=0.5)
        export_parameter_case(current, targets, prefix, history, provenance)
    catch exception
        open(prefix * "-failure.txt", "w") do stream
            showerror(stream, exception, catch_backtrace())
        end
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    parameter_energy_main()
end
