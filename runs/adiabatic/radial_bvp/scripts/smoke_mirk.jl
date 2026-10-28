using BoundaryValueDiffEqMIRK, SciMLBase, NonlinearSolve, ADTypes, Test
println("BVP packages loaded"); flush(stdout)
f! = (du,u,p,t) -> (du[1]=u[2]; du[2]=-p[1]*u[1]; nothing)
left! = (res,u,p) -> (res[1]=u[1]; res[2]=u[2]-1; nothing)
right! = (res,u,p) -> (res[1]=u[1]; nothing)
guess = (p,t) -> [sin(pi*t)/pi,cos(pi*t)]
problem = TwoPointBVProblem(f!,(left!,right!),guess,(0.0,1.0),[9.0];
    bcresid_prototype=(zeros(2),zeros(1)),tune_parameters=true)
@time solution = solve(problem,MIRK4();dt=0.1,abstol=1e-7,reltol=1e-7)
@show solution.retcode solution.prob.p length(solution.t)
@test SciMLBase.successful_retcode(solution)
@test abs(solution.prob.p[1]-pi^2) < 1e-5
@test maximum(abs(solution(t)[1]-sin(pi*t)/pi) for t in range(0,1;length=51)) < 1e-6
