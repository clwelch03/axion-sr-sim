# Rematch the saved populated benchmark on the requested domain, at its
# original per-cloud energies. Run cosine, then chiral from that corrected seed.
include("parameter_energy_target.jl")

function comparison_problem(problem; inner=-80.0, outer=180.0, potential=:cosine)
    isfinite(inner) && isfinite(outer) && inner < problem.integration.match < outer ||
        error("Require ordered finite endpoints around the matching point")
    potential in (:cosine, :chiral) || error("Expected cosine or chiral")
    fields = (; (name => getfield(problem, name) for name in fieldnames(typeof(problem)))...)
    integration = merge(problem.integration, (; inner=Float64(inner), outer=Float64(outer)))
    TwoPhaseProblem(; merge(fields, (; integration, potential))...)
end

function remake_profile_main()
    length(ARGS) in (2, 3) || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY [--check]")
    inputs, output = abspath.(ARGS[1:2])
    checking = length(ARGS) == 3
    checking && ARGS[3] != "--check" && error("Unknown option")
    seed_path = joinpath(inputs, "equal-energy.jls")
    current = deserialize(seed_path)
    @assert current.assessment.passed && current.problem.potential === :cosine
    @assert current.problem.alpha == 0.3 && current.problem.spin == 0.99
    @assert length(current.problem.source_basis.mode_labels) == 14
    config = TOML.parsefile(joinpath(inputs, "request.toml"))
    targets = fill(config["reference_energy"], 2)
    measured = measured_cloud_energies(current; step=0.5)
    @assert maximum(abs.(measured ./ current.energies .- 1)) < 1e-5
    for potential in (:cosine, :chiral)
        problem = comparison_problem(current.problem; inner=config["inner"],
            outer=config["outer"], potential)
        @assert first(current.assessment.positions) <= problem.integration.inner
        @assert last(current.assessment.positions) >= problem.integration.outer
        prediction = energy_predictor(current, current.amplitudes; problem)
        @assert two_phase_profile_admissible(pack_two_phase_profile(prediction.profile, problem), problem)
        @assert all(isfinite, prediction.guess(problem.integration.inner))
        @assert all(isfinite, prediction.guess(problem.integration.outer))
        checking && continue
        mkpath(output)
        prefix = joinpath(output, "joint_" * string(potential))
        any(startswith(name, basename(prefix) * suffix) for name in readdir(output)
            for suffix in (".", "-")) && error("Preserve existing outputs: $prefix")
        history = Dict{String,Any}[]
        provenance = Dict("seed_path"=>seed_path,
            "seed_sha256"=>bytes2hex(sha256(read(seed_path))),
            "request_sha256"=>bytes2hex(sha256(read(joinpath(inputs,"request.toml")))),
            "source_snapshot"=>repository,"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
            "created_utc"=>string(now(UTC)))
        println("BEGIN potential=",potential," domain=",(problem.integration.inner,problem.integration.outer),
            " target energies=",targets); flush(stdout)
        try
            current = solve_parameter_step(current, problem, current.amplitudes, prefix, history)
            current = match_energy_pair(current, targets, prefix, history; final=true)
            current = solve_parameter_step(current, problem, current.amplitudes, prefix, history; tolerance=1e-6)
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
    checking && println("SETUP CHECK PASSED: both potentials, unchanged 14-mode support, original energies")
end

abspath(PROGRAM_FILE) == @__FILE__ && remake_profile_main()
