include("inputs.jl")

function finite_domain_observables(record; stride=2)
    assessment = record.assessment
    case = initial_case(record.case)
    problem = case.problem
    indices = 1:stride:length(assessment.positions)
    branch = (;positions=assessment.positions[indices],
        states=assessment.states[:,indices],spacing=0.25stride)
    result = integrate_two_phase_branch_observables(branch,
        problem.source_basis.mode_labels,problem.spin,problem.alpha,
        real.(assessment.profile.frequencies),problem.cloud_m,
        problem.polar_values,problem.grid,problem.phase1_nodes,problem.phase2_nodes;
        potential=problem.potential)
    return (;energy=result.energy,angular_momentum=result.angular_momentum,
        horizon_energy_flux=-result.lower.radial_energy_flux,
        outer_energy_flux=result.upper.radial_energy_flux,
        horizon_angular_momentum_flux=-result.lower.radial_angular_momentum_flux,
        outer_angular_momentum_flux=result.upper.radial_angular_momentum_flux)
end

function short_reintegration(record)
    case = initial_case(record.case)
    p = case.problem
    a = record.assessment
    rhs! = (du,u,t) -> nonlinear_two_phase_rhs!(du,u,t,p.spin,p.alpha,
        p.components,a.profile.frequencies,p.cloud_m,p.couplings,p.source_basis;
        potential=p.potential)
    defects = zeros(28)
    step_changes = zeros(28)
    for start in -90.0:25.0:285.0
        i = searchsortedfirst(a.positions,start)
        j = searchsortedfirst(a.positions,start+1)
        coarse = integrate_radial(rhs!,a.states[:,i],start,start+1;step=0.05)
        fine = integrate_radial(rhs!,a.states[:,i],start,start+1;step=0.025)
        defects .= max.(defects,abs.(fine-a.states[:,j])./record.state_scales[1:28])
        step_changes .= max.(step_changes,abs.(fine-coarse)./record.state_scales[1:28])
    end
    return (;defects,step_changes)
end

function compare_records(paths)
    records = deserialize.(paths)
    reference = last(records)
    for (path,record) in zip(paths,records)
        println("FILE ",path)
        println("status=",record.retcode," seconds=",record.elapsed,
            " points=",length(record.mesh)," calls=",record.rhs_calls)
        hasproperty(record,:assessment) || continue
        a = record.assessment
        println("frequencies=",a.profile.frequencies,
            " boundary_residual=",a.boundary_residual,
            " defect=",maximum(a.defect_by_state))
        if hasproperty(reference,:assessment) && record.amplitudes == reference.amplitudes
            @assert a.positions == reference.assessment.positions
            denominator = vec(maximum(abs.(reference.assessment.states);dims=2))
            differences = vec(maximum(abs.(a.states-reference.assessment.states);dims=2))
            println("relative_profile_changes_vs_last=",differences./max.(denominator,1e-30))
            println("frequency_changes_vs_last=",a.profile.frequencies-reference.assessment.profile.frequencies)
        end
        quick = get(ENV,"RADIAL_BVP_QUICK_CHECK","false") == "true"
        obs = finite_domain_observables(record;stride=quick ? 4 : 2)
        println(quick ? "finite_domain_observables_h1=" : "finite_domain_observables_h0.5=",obs)
        if !quick
            coarse = finite_domain_observables(record;stride=4)
            println("quadrature_relative_changes_h1_to_h0.5=",
                [(coarse[k]-obs[k])/max(abs(obs[k]),1e-30) for k in keys(obs)])
        end
        println("short_reintegration=",short_reintegration(record))
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_records(ARGS)
end
