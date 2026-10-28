using Test
include("solve_case.jl")
const LS = parentmodule(BackTracking)

@testset "Directional line search on physical collocation residual" begin
    case = initial_case("target_driven")
    setup = augmented_problem(case)
    alg = MIRK4(jac_alg=BVPJacobianAlgorithm(mirk_sparse_ad(28)))
    cache = init(setup.problem,alg;dt=200.0)
    prob = BoundaryValueDiffEqMIRK.__construct_problem(cache,cache.y₀_flat,cache.y₀)
    x = copy(prob.u0)
    fx = similar(x)
    prob.f(fx,x,prob.p)
    J = copy(prob.f.jac_prototype)
    prob.f.jac(J,x,prob.p)
    direct,_,_ = LS.construct_jvp_or_vjp_operator(prob,fx,x;
        autodiff=AutoForwardDiff(;chunksize=1))
    materialized,_,_ = LS.construct_jvp_or_vjp_operator(prob,fx,x;autodiff=nothing)
    for direction in (sin.(eachindex(x)),cos.(eachindex(x)))
        expected = J*direction
        result = similar(fx)
        direct(result,direction,x,prob.p)
        error = norm(result-expected,Inf)/norm(expected,Inf)
        println("Physical JVP relative discrepancy: ",error)
        @test error < 1e-12
        full = similar(fx)
        materialized(full,direction,x,prob.p)
        @test norm(full-result,Inf)/norm(full,Inf) < 1e-12
    end
    println("Operator storage: direct=",Base.summarysize(direct),
        " materialized=",Base.summarysize(materialized))
end

@testset "Backtracking keeps the same acceptance rule" begin
    # Newton's first full step from 0.1 overshoots; both paths must backtrack.
    f!(r,x,p) = (r[1]=x[1]^2-1; nothing)
    jac!(J,x,p) = (J[1,1]=2x[1]; nothing)
    x = [0.1]; fx = [-0.99]; direction = [4.95]
    prob = NonlinearProblem(NonlinearFunction(f!;jac=jac!),x)
    old = init(prob,BackTracking(;maxiters=20),fx,x)
    new = init(prob,BackTracking(;maxiters=20,autodiff=AutoForwardDiff(;chunksize=1)),fx,x)
    a = solve!(old,x,direction)
    b = solve!(new,x,direction)
    @test a.retcode == b.retcode == SciMLBase.ReturnCode.Success
    @test a.step_size ≈ b.step_size
    @test 0 < b.step_size < 1
    println("Accepted step: ",a.step_size," / ",b.step_size)
end
