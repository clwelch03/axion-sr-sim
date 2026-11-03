# Match the Omiya per-cloud energy convention by changing boundary amplitudes
# and solving the full nonlinear equations after every change.
include("collocation_target.jl")
include("comparison_secondary_energy.jl")

function cloud_energy_tables(problem)
    labels = problem.source_basis.mode_labels
    primary_indices = findall(label -> label[2][2] == 0, labels)
    secondary_indices = findall(label -> label[2] == (0,1), labels)
    primary_modes = [(labels[i][1], labels[i][2][1]) for i in primary_indices]
    secondary_modes = [(labels[i][1], 2) for i in secondary_indices]
    grid = angular_grid(length(problem.grid.cosine_theta_nodes), length(problem.phase1_nodes))
    polar = build_polar_values(maximum(first, labels), grid.cosine_theta_nodes)
    # Both energies reconstruct fields; a fundamental-only primary still needs
    # the secondary's m=2 azimuthal harmonic in this shared table.
    maximum_m = max(maximum(last, primary_modes), maximum(last, secondary_modes))
    azimuthal = build_azimuthal_values(maximum_m, grid.azimuth_nodes)
    return (; primary_indices, secondary_indices, primary_modes, secondary_modes,
        grid, polar, azimuthal, mode_count=length(labels))
end

# E1 is the full single-cloud energy of the pure-primary harmonics.
# E2 is quadratic in the secondary on the phase-averaged primary background.
# This is the SAME convention as reference_profiles.jl; it is a comparison
# coordinate, not a unique decomposition of the interacting field's total E.
function comparison_cloud_energies(sampler, frequencies, problem, tables; step=1.0)
    positions = collect(problem.integration.inner:step:problem.integration.outer)
    @assert last(positions) == problem.integration.outer
    primary_rows = vcat(tables.primary_indices, tables.mode_count .+ tables.primary_indices)
    secondary_rows = vcat(tables.secondary_indices, tables.mode_count .+ tables.secondary_indices)
    densities = map(positions) do position
        state = sampler(position)
        primary = state[primary_rows]
        secondary = state[secondary_rows]
        primary_values = Dict(mode => primary[i] for (i,mode) in enumerate(tables.primary_modes))
        e1 = energy_per_tortoise_radius(primary, position, frequencies[1],
            problem.spin, problem.alpha, tables.primary_modes, tables.polar,
            tables.azimuthal, tables.grid; potential=problem.potential)
        e2 = comparison_secondary_energy_density(secondary, position, primary_values,
            frequencies[2], problem.spin, problem.alpha, tables.secondary_modes,
            tables.polar, tables.azimuthal, tables.grid; potential=problem.potential)
        [e1,e2]
    end
    return [integrate_uniform(getindex.(densities,i),step) for i in 1:2]
end

function energy_record(result, assessment, problem, elapsed)
    return (; problem, result.amplitudes, result.state_scales, result.parameter_scales,
        result.packed_reference, result.settings, result.profile, result.successful,
        retcode=string(result.solution.retcode), elapsed,
        mesh=copy(result.solution.t), mesh_states=hcat(result.solution.u...), assessment)
end

function equal_energy_main()
    length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_PREFIX")
    inputs, prefix = abspath(ARGS[1]), abspath(ARGS[2])
    any(ispath(prefix*s) for s in (".jls", ".toml", ".csv", "-history.toml", "-failure.txt")) &&
        error("Preserve existing results; choose a new prefix")
    current = deserialize(joinpath(inputs,"initial.jls"))
    refs = deserialize(joinpath(inputs,"reference-curves.jls"))
    refmeta = TOML.parsefile(joinpath(inputs,"reference-curves.toml"))
    input = deserialize(joinpath(inputs,"joint-input.jls"))
    problem = current.problem
    @assert current.assessment.passed && current.amplitudes == input.cloud_amplitudes
    @assert refmeta["normalization"] == "equal_cloud_energies"
    @assert refmeta["omiya_outer_amplitudes"] == current.amplitudes
    @assert refmeta["input_sha256"]["direct_fig7_20261026_v5-input.jls"] ==
        bytes2hex(sha256(read(joinpath(inputs,"joint-input.jls"))))
    target = refmeta["reference_energy"]
    tables = cloud_energy_tables(problem)
    @assert refs.primary_modes == tables.primary_modes
    @assert refs.secondary_modes == tables.secondary_modes

    # Validate this mode extraction and quadrature against both saved Omiya energies.
    primary_sampler = reference_sampler(refs.omiya_primary_profiles)
    secondary_sampler = reference_sampler(refs.omiya_secondary_profiles)
    function reference_state(position)
        state = zeros(ComplexF64, 2tables.mode_count)
        state[vcat(tables.primary_indices,tables.mode_count .+ tables.primary_indices)] .= primary_sampler(position)
        state[vcat(tables.secondary_indices,tables.mode_count .+ tables.secondary_indices)] .= secondary_sampler(position)
        state
    end
    reference_frequencies = input.profile_guess.frequencies
    reference_energies = comparison_cloud_energies(reference_state,
        reference_frequencies,problem,tables;step=0.5)
    @assert maximum(abs.(reference_energies ./ target .- 1)) < 1e-5
    println("Reference energy check: ", reference_energies, " target=",target); flush(stdout)

    history = Dict{String,Any}[]
    jacobian = 2.0 .* Matrix{Float64}(I,2,2)
    previous_x = previous_y = nothing
    refined = false
    for iteration in 0:8
        assessment = current.assessment
        sampler = position -> hermite_state(assessment.positions,assessment.states,position)
        energies = comparison_cloud_energies(sampler,assessment.profile.frequencies,problem,tables)
        error_vector = log.(energies ./ target)
        x = log.(current.amplitudes)
        push!(history,Dict("iteration"=>iteration,"amplitudes"=>current.amplitudes,
            "energies"=>energies,"relative_energy_error"=>energies ./ target .- 1,
            "differential_tolerance"=>current.settings.tolerance,
            "boundary_residual"=>assessment.boundary_residual,
            "defect"=>maximum(assessment.defect_by_state)))
        open(prefix*"-history.toml","w") do stream
            TOML.print(stream,Dict("target_energy"=>target,"steps"=>history);sorted=true)
        end
        println("ENERGY iteration=",iteration," A=",current.amplitudes,
            " E=",energies," relative_error=",energies ./ target .- 1); flush(stdout)

        if maximum(abs,error_vector) < 3e-5
            if refined
                fine_energies = comparison_cloud_energies(sampler,
                    assessment.profile.frequencies,problem,tables;step=0.5)
                @assert maximum(abs.(fine_energies ./ energies .- 1)) < 5e-6
                @assert maximum(abs.(fine_energies ./ target .- 1)) < 1e-4
                current = merge(current,(; energies=fine_energies, reference_energy=target))
                break
            end
            next_amplitudes = current.amplitudes
            refined = true
        else
            if previous_x !== nothing
                dx, dy = x-previous_x, error_vector-previous_y
                dot(dx,dx) > 1e-20 && (jacobian += (dy-jacobian*dx)*transpose(dx)/dot(dx,dx))
            end
            correction = -(jacobian \ error_vector)
            correction *= min(1.0,0.2/maximum(abs,correction))
            next_amplitudes = exp.(x+correction)
            previous_x, previous_y = x,error_vector
        end
        iteration < 8 || error("Energy targeting exceeded eight solves")
        ratios = next_amplitudes ./ current.amplitudes
        factors = [prod(ratios .^ abs.(collect(q))) for (l,q) in problem.source_basis.mode_labels]
        state_factors = vcat(factors,factors)
        guess = position -> sampler(position) .* state_factors
        profile = TwoPhaseProfile(copy(assessment.profile.frequencies),
            assessment.profile.horizon_coefficients .* factors,
            assessment.profile.outer_coefficients .* factors)
        scales = two_phase_collocation_scales(problem,next_amplitudes,profile,guess)
        tolerance = refined ? 1e-6 : 1e-5
        timing = @timed solve_two_phase_collocation(problem,next_amplitudes,profile,guess;
            scales...,initial_step=refined ? 0.5 : 1.0,tolerance,
            nonlinear_tolerance=tolerance/100,ad_chunksize=28,time_limit=1800.0,trace=true)
        result = timing.value
        assessment = assess_two_phase_collocation(result,problem)
        current = energy_record(result,assessment,problem,timing.time)
        serialize(prefix*"-step$(iteration+1).jls",current)
        assessment.passed || error("Radial assessment failed at energy iteration $(iteration+1)")
    end
    hasproperty(current,:energies) || error("No accepted energy-matched solution")
    assessment = current.assessment
    # The total joint energy includes every retained mode and the full potential.
    indices = 1:4:length(assessment.positions)
    branch = (;positions=assessment.positions[indices],states=assessment.states[:,indices],spacing=1.0)
    total = integrate_two_phase_branch_observables(branch,problem.source_basis.mode_labels,
        problem.spin,problem.alpha,real.(assessment.profile.frequencies),problem.cloud_m,
        problem.polar_values,problem.grid,problem.phase1_nodes,problem.phase2_nodes;potential=problem.potential)
    serialize(prefix*".jls",merge(current,(;history,total_observables=total)))
    metadata = Dict("passed"=>true,"normalization"=>"equal_cloud_energies",
        "energy_convention"=>"Omiya: isolated cosine primary; quadratic secondary on averaged primary; mixed modes excluded from per-cloud labels",
        "reference_energy"=>target,"joint_cloud_energies"=>current.energies,
        "joint_total_energy"=>total.energy,"joint_total_angular_momentum"=>total.angular_momentum,
        "cloud_amplitudes"=>current.amplitudes,"original_outer_amplitudes"=>input.cloud_amplitudes,
        "spin"=>problem.spin,"alpha"=>problem.alpha,"potential"=>string(problem.potential),
        "inner"=>problem.integration.inner,"outer"=>problem.integration.outer,
        "mode_labels"=>string.(problem.source_basis.mode_labels),
        "polar_nodes"=>length(problem.grid.cosine_theta_nodes),
        "phase_nodes"=>[length(problem.phase1_nodes),length(problem.phase2_nodes)],
        "horizon_order"=>problem.integration.horizon_order,"outer_order"=>problem.integration.outer_order,
        "mesh_points"=>length(current.mesh),"boundary_residual"=>assessment.boundary_residual,
        "parameter_drift"=>assessment.parameter_drift,"independent_off_mesh_defect"=>maximum(assessment.defect_by_state),
        "requested_tolerance"=>current.settings.tolerance,
        "frequencies_real"=>real.(assessment.profile.frequencies),
        "frequencies_imaginary"=>imag.(assessment.profile.frequencies),
        "reference_curve_metadata"=>refmeta,"created_utc"=>string(now(UTC)),
        "source_snapshot"=>repository,"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "input_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(inputs,name)))) for name in readdir(inputs)))
    open(prefix*".toml","w") do stream
        TOML.print(stream,metadata;sorted=true)
    end
    names = (:omiya_primary,:omiya_secondary,:linear_primary,:linear_secondary)
    samplers = [reference_sampler(getproperty(refs,Symbol(name,"_profiles"))) for name in names]
    open(prefix*".csv","w") do stream
        println(stream,"curve,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for (j,position) in enumerate(assessment.positions)
            for (name,sample) in zip(names,samplers)
                state = sample(position)
                R,Q = state[1],state[length(state)÷2+1]
                ell = endswith(string(name),"primary") ? 1 : 2
                @printf(stream,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",name,ell,ell,position,real(R),imag(R),abs(R),real(Q),imag(Q))
            end
            for (name,i,ell) in (("joint_primary",first(tables.primary_indices),1),("joint_secondary",first(tables.secondary_indices),2))
                R,Q = assessment.states[i,j],assessment.states[tables.mode_count+i,j]
                @printf(stream,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",name,ell,ell,position,real(R),imag(R),abs(R),real(Q),imag(Q))
            end
        end
    end
    println("EQUAL ENERGY ACCEPTED: E=",current.energies," total E=",total.energy,
        " amplitudes=",current.amplitudes); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        equal_energy_main()
    catch exception
        if length(ARGS)==2
            open(abspath(ARGS[2])*"-failure.txt","w") do stream
                showerror(stream,exception,catch_backtrace())
            end
        end
        rethrow()
    end
end
