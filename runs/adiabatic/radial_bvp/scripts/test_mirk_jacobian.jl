using BoundaryValueDiffEqMIRK, ForwardDiff, LinearAlgebra, Test
include("mirk_coupling_patch.jl")

function check_jacobian()
    A = reshape(sin.(1.0:36.0),6,6)./10
    exact(t) = exp(A*t)*ones(6)
    f! = (du,u,p,t) -> (mul!(du,A,u);nothing)
    left! = (res,u,p) -> (res .= u[1:3].-1;nothing)
    right! = (res,u,p) -> (res .= u[4:6].-exact(1.0)[4:6];nothing)
    guess = (p,t) -> exact(t).+0.01sin(t)
    prob = TwoPointBVProblem(f!,(left!,right!),guess,(0.,1.);
        bcresid_prototype=(zeros(3),zeros(3)))
    alg = MIRK4(jac_alg=BVPJacobianAlgorithm(
        mirk_sparse_ad(12)))
    cache = init(prob,alg;dt=0.5)
    nlprob = BoundaryValueDiffEqMIRK.__construct_problem(cache,cache.y₀_flat,cache.y₀)
    x = copy(nlprob.u0)
    J = copy(nlprob.f.jac_prototype)
    nlprob.f.jac(J,x,nlprob.p)
    direct = ForwardDiff.jacobian(x) do z
        out = similar(z)
        nlprob.f(out,z,nlprob.p)
        out
    end
    error = norm(Matrix(J)-direct,Inf)
    println("Jacobian max entry error: ",error);flush(stdout)
    return error,prob,alg,exact
end

unpatched,_,_,_ = check_jacobian()
@test unpatched > 1e-4
apply_mirk_coupling_patch!()
patched,prob,alg,exact = check_jacobian()
@test patched < 1e-13
sol = solve(prob,alg;dt=0.2,abstol=1e-7)
@test successful_retcode(sol)
@test maximum(norm(sol(t)-exact(t),Inf) for t in range(0,1;length=51)) < 1e-6
println("Coupled linear BVP and full Jacobian checks passed.")
