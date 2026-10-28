using Profile, TOML
include("solve_case.jl")

mutable struct CallbackTimings
    rhs_calls::Dict{Int,Int}
    rhs_ns::Dict{Int,UInt64}
    bc_calls::Int
    bc_ns::UInt64
end
CallbackTimings() = CallbackTimings(Dict(),Dict(),0,UInt64(0))

function measured_problem(case)
    setup = augmented_problem(case;time_limit=3500.0)
    metrics = CallbackTimings()
    function rhs!(du,u,p,t)
        width = eltype(u) <: ForwardDiff.Dual ? ForwardDiff.npartials(eltype(u)) : 0
        started = time_ns()
        setup.problem.f(du,u,p,t)
        metrics.rhs_ns[width] = get(metrics.rhs_ns,width,UInt64(0))+time_ns()-started
        metrics.rhs_calls[width] = get(metrics.rhs_calls,width,0)+1
        nothing
    end
    function left!(r,u,p)
        started = time_ns()
        setup.cb.bca!(r,u,p)
        metrics.bc_ns += time_ns()-started
        metrics.bc_calls += 1
        nothing
    end
    function right!(r,u,p)
        started = time_ns()
        setup.cb.bcb!(r,u,p)
        metrics.bc_ns += time_ns()-started
        metrics.bc_calls += 1
        nothing
    end
    problem = TwoPointBVProblem(rhs!,(left!,right!),setup.problem.u0,(-100.0,300.0);
        bcresid_prototype=(zeros(56),zeros(56)))
    (;problem,metrics,cb=setup.cb)
end

function speed_algorithm(variant)
    linesearch = if variant == "baseline"
        BackTracking(;maxiters=20)
    elseif variant == "directional"
        BackTracking(;maxiters=20,autodiff=AutoForwardDiff(;chunksize=1))
    else
        error("Unknown speed variant $variant")
    end
    chunksize = parse(Int,get(ENV,"RADIAL_BVP_AD_CHUNKSIZE","28"))
    MIRK4(;jac_alg=BVPJacobianAlgorithm(mirk_sparse_ad(chunksize)),
        nlsolve=NewtonRaphson(;jacobian_reuse=false,linesearch),max_num_subintervals=1200)
end

function benchmark_run(case,variant,label,prefix;profile=false)
    output = prefix*"-"*label
    isfile(output*".jls") && error("Existing benchmark $output")
    setup = measured_problem(case)
    alg = speed_algorithm(variant)
    GC.gc()
    Profile.clear()
    println("BEGIN variant=",variant," pass=",label);flush(stdout)
    timing = if profile
        @timed Profile.@profile solve(setup.problem,alg;dt=1.0,abstol=1e-4,
            nlsolve_kwargs=(;abstol=1e-6,reltol=0.0,maxiters=20,maxtime=1800.,show_trace=Val(true)))
    else
        @timed solve(setup.problem,alg;dt=1.0,abstol=1e-4,
            nlsolve_kwargs=(;abstol=1e-6,reltol=0.0,maxiters=20,maxtime=1800.,show_trace=Val(true)))
    end
    sol = timing.value
    println("END variant=",variant," pass=",label," seconds=",timing.time,
        " compile_seconds=",timing.compile_time," retcode=",sol.retcode);flush(stdout)
    if profile
        data,lineinfo = Profile.retrieve()
        serialize(output*"-profile.jls",(;data,lineinfo))
        open(output*"-profile-flat.txt","w") do io
            Profile.print(IOContext(io,:displaysize=>(1000,240)),data,lineinfo;
                format=:flat,sortedby=:count,mincount=10,C=true,threads=1)
        end
        open(output*"-profile-tree.txt","w") do io
            Profile.print(IOContext(io,:displaysize=>(1000,240)),data,lineinfo;
                format=:tree,mincount=50,C=false,threads=1)
        end
    end
    assessment_seconds = @elapsed assessment = assess_solution(sol,case,setup.cb)
    passed_residual_checks = successful_retcode(sol) && all(isfinite,assessment.states) &&
        assessment.boundary_residual<=1e-6 && maximum(assessment.defect_by_state)<=1e-4
    metrics = setup.metrics
    chunksize = parse(Int,get(ENV,"RADIAL_BVP_AD_CHUNKSIZE","28"))
    report = Dict("variant"=>variant,"pass"=>label,"profiled"=>profile,
        "ad_chunksize"=>chunksize,
        "elapsed_seconds"=>timing.time,"compilation_seconds"=>timing.compile_time,
        "recompilation_seconds"=>timing.recompile_time,"gc_seconds"=>timing.gctime,
        "allocated_bytes"=>timing.bytes,"assessment_seconds"=>assessment_seconds,
        "rhs_calls"=>Dict(string(k)=>v for (k,v) in metrics.rhs_calls),
        "rhs_seconds"=>Dict(string(k)=>v/1e9 for (k,v) in metrics.rhs_ns),
        "boundary_calls"=>metrics.bc_calls,"boundary_seconds"=>metrics.bc_ns/1e9,
        "mesh_points"=>length(sol.t),"retcode"=>string(sol.retcode),
        "hostname"=>gethostname(),"cpu"=>Sys.cpu_info()[1].model,
        "julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),
        "boundary_residual"=>assessment.boundary_residual,
        "off_mesh_defect"=>maximum(assessment.defect_by_state),
        "passed_residual_checks"=>passed_residual_checks)
    open(output*".toml","w") do io
        TOML.print(io,report)
    end
    settings = (;tolerance=1e-4,dt=1.,algorithm="MIRK4",variant,ad_chunksize=chunksize)
    source_hashes=Dict(f=>file_hash(joinpath(@__DIR__,f)) for f in
        ("inputs.jl","solve_case.jl","speed_benchmark.jl","mirk_coupling_patch.jl","mirk_remesh_patch.jl"))
    record = (;case=case.name,settings,case.provenance,source_hashes,case.amplitudes,
        case.state_scales,case.parameter_scales,case.packed_reference,
        retcode=string(sol.retcode),successful=successful_retcode(sol),
        elapsed=timing.time,compilation_seconds=timing.compile_time,
        rhs_calls=[get(metrics.rhs_calls,0,0),sum(v for (k,v) in metrics.rhs_calls if k!=0)],
        mesh=copy(sol.t),mesh_states=hcat(sol.u...),assessment,passed_residual_checks)
    serialize(output*".jls",record)
    passed_residual_checks || error("Benchmark did not pass residual checks")
    println("METRICS ",report);flush(stdout)
    nothing
end

function speed_main(args)
    length(args)==2 || error("Usage: speed_benchmark.jl VARIANT OUTPUT_PREFIX")
    variant,prefix=args
    case=initial_case("target_driven")
    Profile.init(n=20_000_000,delay=0.01)
    for label in ("cold","warm")
        Base.invokelatest(benchmark_run,case,variant,label,prefix;profile=label=="warm")
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    speed_main(ARGS)
end
