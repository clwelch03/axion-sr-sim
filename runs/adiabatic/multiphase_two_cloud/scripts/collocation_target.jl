# Reproduce the refined target with current src and export six comparable curves.
# Usage: julia --project=. collocation_target.jl INPUT_DIRECTORY OUTPUT_PREFIX
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
const repository = normpath(joinpath(@__DIR__, "../../../.."))
include(joinpath(repository, "src/adiabatic/primary_secondary/load.jl"))
include(joinpath(repository, "src/adiabatic/multiphase/load_collocation.jl"))
if !isdefined(@__MODULE__, :assess_two_phase_collocation)
    println("Using the tested assessment block packaged with this driver.")
    include(joinpath(@__DIR__, "assess_collocation_target.jl"))
end

"""
    hermite_state(positions, states, position)

Return physical complex `[R; Q]` at a covered tortoise radius using cubic Hermite
interpolation. Supply increasing sample positions and columns ordered as all
fields followed by all slopes. The returned Q is the derivative of the
interpolated R; input arrays are read without mutation.
"""
function hermite_state(positions, states, position)
    first(positions) <= position <= last(positions) || error("Guess does not cover radius $position")
    count = size(states, 1) ÷ 2
    index = clamp(searchsortedlast(positions, position), 1, length(positions)-1)
    width = positions[index+1] - positions[index]
    fraction = (position - positions[index]) / width
    left, right = states[1:count,index], states[1:count,index+1]
    left_slope, right_slope = states[count+1:end,index], states[count+1:end,index+1]
    fields = (2fraction^3-3fraction^2+1).*left +
        (fraction^3-2fraction^2+fraction)*width.*left_slope +
        (-2fraction^3+3fraction^2).*right + (fraction^3-fraction^2)*width.*right_slope
    slopes = (6fraction^2-6fraction)/width.*left +
        (3fraction^2-4fraction+1).*left_slope +
        (-6fraction^2+6fraction)/width.*right + (3fraction^2-2fraction).*right_slope
    return vcat(fields, slopes)
end

"""
    reference_sampler(radial)

Return a callable physical `[R; Q]` interpolator for two saved shooting branches.
Select the horizon branch at or below the matching radius and the outer branch
above it, preserving their separate values at the join. Orient the copied outer
samples toward increasing radius and use [`hermite_state`](@ref) within a branch.
"""
function reference_sampler(radial)
    horizon = radial.from_horizon
    outer = (; positions=reverse(radial.from_outer.positions), states=reverse(radial.from_outer.states; dims=2))
    return position -> begin
        branch = position <= last(horizon.positions) ? horizon : outer
        hermite_state(branch.positions, branch.states, position)
    end
end

function main()
    length(ARGS) == 2 || error("Usage: INPUT_DIRECTORY OUTPUT_PREFIX")
    inputs, prefix = abspath(ARGS[1]), abspath(ARGS[2])
    any(ispath(prefix * suffix) for suffix in (".jls", ".toml", ".csv", "-setup.jls", "-failure.txt")) &&
        error("Preserve existing output; choose a new prefix")
    mkpath(dirname(prefix))
    input = deserialize(joinpath(inputs, "joint-input.jls"))
    old = deserialize(joinpath(inputs, "initial.jls"))
    references = deserialize(joinpath(inputs, "reference-curves.jls"))
    reference_metadata = TOML.parsefile(joinpath(inputs, "reference-curves.toml"))
    problem, amplitudes = input.problem, input.cloud_amplitudes
    @assert old.successful && old.amplitudes == amplitudes
    @assert old.assessment.boundary_residual <= old.settings.tolerance / 100
    @assert maximum(old.assessment.defect_by_state) <= old.settings.tolerance
    @assert old.assessment.parameter_drift <= 1e-10
    @assert reference_metadata["normalization"] == "fixed_outer_coefficients"
    @assert reference_metadata["cloud_amplitudes"] == amplitudes
    @assert reference_metadata["spin"] == problem.spin && reference_metadata["alpha"] == problem.alpha
    @assert reference_metadata["input_sha256"]["direct_fig7_20261026_v5-input.jls"] ==
        bytes2hex(sha256(read(joinpath(inputs, "joint-input.jls"))))
    mode_count = length(problem.source_basis.mode_labels)
    @assert mode_count == 14 && size(old.assessment.states,1) == 2mode_count
    initial_guess = position -> hermite_state(old.assessment.positions, old.assessment.states, position)
    reference_profile = old.assessment.profile
    scales = two_phase_collocation_scales(problem, amplitudes, reference_profile, initial_guess)
    serialize(prefix * "-setup.jls", (; problem, amplitudes, reference_profile, scales,
        input_hashes=Dict(name => bytes2hex(sha256(read(joinpath(inputs,name))))
            for name in readdir(inputs))))
    println("Current-src target: modes=",mode_count," amplitudes=",amplitudes,
        " tolerance=1e-6 initial_step=0.5 chunksize=28"); flush(stdout)
    timing = @timed solve_two_phase_collocation(problem, amplitudes, reference_profile, initial_guess;
        scales..., initial_step=0.5, tolerance=1e-6, nonlinear_tolerance=1e-8,
        ad_chunksize=28, time_limit=2400.0, trace=true)
    result = timing.value
    record = (; problem, amplitudes, result.state_scales, result.parameter_scales,
        result.packed_reference, result.settings, result.profile, result.successful,
        retcode=string(result.solution.retcode), elapsed=timing.time,
        compilation_seconds=timing.compile_time, mesh=copy(result.solution.t),
        mesh_states=hcat(result.solution.u...), slurm_job_id=get(ENV,"SLURM_JOB_ID","local"))
    serialize(prefix * ".jls.tmp", record)
    mv(prefix * ".jls.tmp", prefix * ".jls")
    assessment = assess_two_phase_collocation(result, problem)
    frequency_change = maximum(abs.(assessment.profile.frequencies - reference_profile.frequencies))
    @assert assessment.positions == old.assessment.positions
    relative_profile_changes = vec(maximum(abs.(assessment.states-old.assessment.states); dims=2)) ./
        max.(vec(maximum(abs.(old.assessment.states); dims=2)), 1e-30)
    record = merge(record, (; assessment, frequency_change, relative_profile_changes))
    serialize(prefix * ".jls.tmp", record)
    mv(prefix * ".jls.tmp", prefix * ".jls"; force=true)
    metadata = Dict("passed"=>assessment.passed, "normalization"=>"fixed_outer_coefficients",
        "cloud_amplitudes"=>amplitudes, "spin"=>problem.spin, "alpha"=>problem.alpha,
        "potential"=>string(problem.potential), "mode_labels"=>string.(problem.source_basis.mode_labels),
        "polar_nodes"=>length(problem.grid.cosine_theta_nodes),
        "phase_nodes"=>[length(problem.phase1_nodes),length(problem.phase2_nodes)],
        "inner"=>problem.integration.inner,"outer"=>problem.integration.outer,
        "horizon_order"=>problem.integration.horizon_order,"outer_order"=>problem.integration.outer_order,
        "mesh_points"=>length(result.solution.t), "retcode"=>record.retcode,
        "elapsed_seconds"=>timing.time,"compilation_seconds"=>timing.compile_time,
        "boundary_residual"=>assessment.boundary_residual,"parameter_drift"=>assessment.parameter_drift,
        "independent_off_mesh_defect"=>maximum(assessment.defect_by_state),
        "requested_tolerance"=>result.settings.tolerance,"frequency_change_vs_refined"=>frequency_change,
        "relative_profile_changes_vs_refined"=>relative_profile_changes,
        "frequencies_real"=>real.(assessment.profile.frequencies),
        "frequencies_imaginary"=>imag.(assessment.profile.frequencies),
        "slurm_job_id"=>record.slurm_job_id,"created_utc"=>string(now(UTC)),
        "source_snapshot"=>repository,"reference_curve_metadata"=>reference_metadata)
    open(prefix * ".toml", "w") do stream
        TOML.print(stream, metadata; sorted=true)
    end
    assessment.passed || error("Physical assessment failed; inspect saved diagnostics")
    frequency_change < 1e-7 || error("Target frequencies changed unexpectedly; inspect branch consistency")

    # Export every curve on the same grid, preserving the reference's join.
    reference_names = (:omiya_primary, :omiya_secondary, :linear_primary, :linear_secondary)
    samplers = [reference_sampler(getproperty(references, Symbol(name,"_profiles"))) for name in reference_names]
    labels = problem.source_basis.mode_labels
    primary_index = findfirst(==((1,(1,0))), labels)
    secondary_index = findfirst(==((2,(0,1))), labels)
    @assert primary_index !== nothing && secondary_index !== nothing
    open(prefix * ".csv.tmp", "w") do stream
        println(stream,"curve,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for (index, position) in enumerate(assessment.positions)
            for (name, sampler) in zip(reference_names, samplers)
                state = sampler(position)
                count = length(state) ÷ 2
                field, slope = state[1], state[count+1]
                ell = endswith(string(name),"primary") ? 1 : 2
                @printf(stream,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    string(name),ell,ell,position,real(field),imag(field),abs(field),real(slope),imag(slope))
            end
            for (name, mode_index, ell) in (("joint_primary",primary_index,1),("joint_secondary",secondary_index,2))
                field = assessment.states[mode_index,index]
                slope = assessment.states[mode_count+mode_index,index]
                @printf(stream,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    name,ell,ell,position,real(field),imag(field),abs(field),real(slope),imag(slope))
            end
        end
    end
    mv(prefix * ".csv.tmp", prefix * ".csv")
    println("TARGET ACCEPTED: boundary=",assessment.boundary_residual,
        " defect=",maximum(assessment.defect_by_state)," frequency_change=",frequency_change)
    println("Exported six physical curves: ",prefix,".csv"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch exception
        if length(ARGS) == 2
            open(abspath(ARGS[2]) * "-failure.txt", "w") do stream
                showerror(stream,exception,catch_backtrace())
            end
        end
        rethrow()
    end
end
