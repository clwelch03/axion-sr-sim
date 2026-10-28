using BoundaryValueDiffEqMIRK,Test
include("mirk_remesh_patch.jl")

function remesh_error()
    f! = (du,u,p,t) -> (du .= 1;nothing)
    left! = (res,u,p) -> (res[1]=u[1];nothing)
    right! = (res,u,p) -> (res[1]=u[2]-1;nothing)
    problem = TwoPointBVProblem(f!,(left!,right!),(p,t)->[t,t],(0.,1.);
        bcresid_prototype=(zeros(1),zeros(1)))
    cache = init(problem,MIRK4();dt=0.5)
    # Exact stages for y'=1, y(0)=0, so the interpolant must be y=t.
    for stages in cache.k_discrete
        fill!(hasproperty(stages,:du) ? stages.du : stages,1.0)
    end
    for stages in cache.k_interp.u
        fill!(stages,1.0)
    end
    old_mesh,old_dt = copy(cache.mesh),copy(cache.mesh_dt)
    new_mesh = collect(0.0:0.1:1.0)
    append!(cache.y₀.u,[zeros(2) for _ in 1:length(new_mesh)-length(cache.y₀.u)])
    for (i,t) in enumerate(new_mesh)
        BoundaryValueDiffEqMIRK.interp_eval!(cache.y₀.u[i],cache,t,old_mesh,old_dt)
    end
    maximum(abs(cache.y₀.u[i][1]-t) for (i,t) in enumerate(new_mesh))
end

before = remesh_error()
apply_mirk_remesh_patch!()
after = remesh_error()
println("Exact y=t remesh error before correction: ",before)
println("Exact y=t remesh error after correction: ",after)
@test before > 0.1
@test after < 1e-14
