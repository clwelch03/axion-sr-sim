using Test
include("solve_case.jl")
case = initial_case("weak")
setup = augmented_problem(case)
algorithm = MIRK4(jac_alg=BVPJacobianAlgorithm(
    mirk_sparse_ad(28)))
# Three nodes suffice to exercise both boundaries and every block coupling.
# This is a derivative check, not an attempted physical discretization.
cache = init(setup.problem,algorithm;dt=200.0)
nlprob = BoundaryValueDiffEqMIRK.__construct_problem(cache,cache.y₀_flat,cache.y₀)
x = copy(nlprob.u0)
J = copy(nlprob.f.jac_prototype)
nlprob.f.jac(J,x,nlprob.p)
direct = ForwardDiff.jacobian(x) do z
    out = similar(z)
    nlprob.f(out,z,nlprob.p)
    out
end
discrepancy = norm(Matrix(J)-direct,Inf)/norm(direct,Inf)
println("Physical collocation Jacobian vs uncompressed AD: ",discrepancy)
@test discrepancy < 1e-12
