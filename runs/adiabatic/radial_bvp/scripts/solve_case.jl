using BoundaryValueDiffEqMIRK, SciMLBase, NonlinearSolve, ADTypes
include("inputs.jl")
include("mirk_coupling_patch.jl")
apply_mirk_coupling_patch!()
include("mirk_remesh_patch.jl")
apply_mirk_remesh_patch!()

# Explicit augmentation avoids relying on the package's parameter/BC wrapper.
# The corrected two-point Jacobian band retains all local mode couplings.
# No boundary formula or radial source is changed.
function augmented_problem(case; restart=nothing, time_limit=2400.0)
    cb = collocation_callbacks(case)
    started = time()
    heartbeat = Ref(started)
    counts = [0,0]
    function rhs!(du,u,p,t)
        counts[eltype(u) <: ForwardDiff.Dual ? 2 : 1] += 1
        now = time()
        now-started <= time_limit || error("Experiment wall-time limit reached")
        if now-heartbeat[] > 30
            println("heartbeat seconds=",round(now-started;digits=1)," rhs_calls=",counts)
            flush(stdout)
            heartbeat[] = now
        end
        cb.rhs!(du,u,view(u,57:112),t)
        du[57:112] .= 0
        nothing
    end
    guess = (p,t) -> vcat(cb.guess(zeros(56),t),zeros(56))
    if restart !== nothing
        old = deserialize(restart)
        old.amplitudes == case.amplitudes || error("Restart amplitudes differ")
        old.state_scales == case.state_scales || error("Restart scales differ")
        old.parameter_scales == case.parameter_scales || error("Restart parameter scales differ")
        old.packed_reference == case.packed_reference || error("Restart parameter reference differs")
        guess = (p,t) -> begin
            i = clamp(searchsortedlast(old.mesh,t),1,length(old.mesh)-1)
            w = (t-old.mesh[i])/(old.mesh[i+1]-old.mesh[i])
            (1-w).*old.mesh_states[:,i] .+ w.*old.mesh_states[:,i+1]
        end
    end
    problem = TwoPointBVProblem(rhs!,(cb.bca!,cb.bcb!),guess,(-100.0,300.0);
        bcresid_prototype=(zeros(56),zeros(56)))
    return (; problem, cb, counts, started)
end

function assess_solution(sol,case,cb)
    p = sol.u[1][57:112]
    parameter_drift = maximum(norm(u[57:112]-p,Inf) for u in sol.u)
    boundary = zeros(112)
    cb.bca!(view(boundary,1:56),first(sol.u),nothing)
    cb.bcb!(view(boundary,57:112),last(sol.u),nothing)
    defect_by_state = zeros(56)
    absolute_defect_by_state = zeros(56)
    # These points are distinct from MIRK's defect-control sample points.
    for i in 1:length(sol.t)-1, fraction in (0.211324865405187,0.788675134594813)
        t = sol.t[i]+fraction*(sol.t[i+1]-sol.t[i])
        u,du = sol(t),sol(t,Val{1})
        f = zeros(56)
        cb.rhs!(f,u,view(u,57:112),t)
        defect_by_state .= max.(defect_by_state,abs.(du[1:56]-f)./(1 .+ abs.(f)))
        absolute_defect_by_state .= max.(absolute_defect_by_state,
            case.state_scales.*abs.(du[1:56]-f))
    end
    positions = collect(-100.0:0.25:300.0)
    states = hcat([unpack_complex(case.state_scales.*sol(t)[1:56]) for t in positions]...)
    profile = cb.profile_at(p)
    reference = hcat([case.guess(t) for t in positions]...)
    profile_change_by_mode = vec(maximum(abs.(states-reference);dims=2)) ./
        case.state_scales[1:28]
    return (; profile, positions, states, parameter_drift,
        boundary_residual=norm(boundary,Inf), defect_by_state,
        absolute_defect_by_state, profile_change_by_mode)
end

function main(args;case_override=nothing)
    length(args) >= 4 || error("Usage: solve_case.jl CASE TOL INITIAL_DT OUTPUT_PREFIX [RESTART]")
    name, tolerance, dt, prefix = args[1],parse(Float64,args[2]),parse(Float64,args[3]),args[4]
    any(isfile(prefix*suffix) for suffix in (".jls","-failure.txt","-setup.jls")) &&
        error("Refusing to overwrite existing experiment $prefix")
    mkpath(dirname(prefix))
    case = case_override === nothing ? initial_case(name) : case_override
    use_linesearch = get(ENV,"RADIAL_BVP_LINESEARCH","true") == "true"
    settings = (; tolerance,dt,algorithm="MIRK4",
        nonlinear_algorithm=use_linesearch ? "NewtonRaphson+BackTracking" : "NewtonRaphson",
        ad_chunksize=28,coloring="GreedyColoringAlgorithm",
        max_num_subintervals=1200,nonlinear_maxiters=20,time_limit=2400.0)
    source_hashes = Dict(f=>file_hash(joinpath(@__DIR__,f)) for f in
        ("inputs.jl","solve_case.jl","mirk_coupling_patch.jl","mirk_remesh_patch.jl"))
    serialize(prefix*"-setup.jls",(; settings,case=case.name,case.provenance,
        case.amplitudes,case.state_scales,case.parameter_scales,case.packed_reference,source_hashes))
    println("case=",name," settings=",settings," amplitudes=",case.amplitudes)
    flush(stdout)
    setup = augmented_problem(case;restart=length(args)>4 ? args[5] : nothing,
        time_limit=settings.time_limit)
    jac = BVPJacobianAlgorithm(mirk_sparse_ad(28))
    nonlinear_solver = use_linesearch ?
        NewtonRaphson(;jacobian_reuse=false,linesearch=BackTracking(;maxiters=20)) :
        NewtonRaphson(;jacobian_reuse=false)
    algorithm = MIRK4(;jac_alg=jac,nlsolve=nonlinear_solver,
        max_num_subintervals=settings.max_num_subintervals)
    timing = @timed solve(setup.problem,algorithm;dt,abstol=tolerance,
        nlsolve_kwargs=(;abstol=tolerance/100,reltol=0.0,maxiters=settings.nonlinear_maxiters,
            maxtime=900.0,show_trace=Val(true)))
    sol,elapsed = timing.value,timing.time
    println("retcode=",sol.retcode," seconds=",elapsed,
        " compilation_seconds=",timing.compile_time," mesh_points=",length(sol.t))
    flush(stdout)
    # Save even unsuccessful returns, explicitly labelled by return code.
    preliminary = (;case=name,settings,case.provenance,source_hashes,
        case.amplitudes,case.state_scales,case.parameter_scales,case.packed_reference,
        retcode=string(sol.retcode),successful=SciMLBase.successful_retcode(sol),
        elapsed,compilation_seconds=timing.compile_time,recompilation_seconds=timing.recompile_time,
        rhs_calls=copy(setup.counts),mesh=copy(sol.t),mesh_states=hcat(sol.u...))
    serialize(prefix*".jls",preliminary)
    assessment = assess_solution(sol,case,setup.cb)
    passed_residual_checks = preliminary.successful && all(isfinite,assessment.states) &&
        assessment.boundary_residual <= tolerance/100 &&
        maximum(assessment.defect_by_state) <= tolerance
    record = merge(preliminary,(;assessment,passed_residual_checks))
    serialize(prefix*".jls",record)
    println("frequencies=",assessment.profile.frequencies,
        " boundary_residual=",assessment.boundary_residual,
        " off_mesh_defect=",maximum(assessment.defect_by_state),
        " parameter_drift=",assessment.parameter_drift)
    println("profile_change_by_mode=",assessment.profile_change_by_mode)
    flush(stdout)
    return record
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch exception
        if length(ARGS) >= 4 && !isfile(ARGS[4]*"-failure.txt")
            open(ARGS[4]*"-failure.txt","w") do io
                showerror(io,exception,catch_backtrace())
            end
        end
        rethrow()
    end
end
