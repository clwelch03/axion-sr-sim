using Test
include("run_case.jl")

root = normpath(joinpath(@__DIR__, "../../.."))
directory = joinpath(root, "runs/adiabatic/radialsolvertest")
saved = deserialize(joinpath(directory, "seed_cosine.jls"))
context = benchmark_context(saved)
direct = deserialize(joinpath(directory, "shooting_ad.jls"))
deuflhard = deserialize(joinpath(directory, "shooting_ad_deuflhard.jls"))
results = NamedTuple[]

@testset "AD shooting: subdivision and correction-strategy agreement" begin
    @test direct.result.status == deuflhard.result.status == "converged"
    history = copy(saved.continuation)
    for target in (last(history).amplitude+0.05, direct.amplitude)
        prediction = predict_nonlinear_mode_cubic(history, target; history_spacing=0.025)
        problem = shooting_problem(context, target, prediction)
        jacobian = state -> ForwardDiff.jacobian(problem.residual, state)
        result = correct_system(problem.residual, jacobian, problem.initial;
            admissible=state -> all(isfinite,state) && 0 < state[1] < context.alpha)
        @test result.status == "converged"
        solution = physical_solution(result.state, target, length(result.history)-1, result.residual)
        diagnostics = evaluate_solution(solution, target, context)
        @test diagnostics.replay_residual <= 1e-9
        push!(history, (; amplitude=target, solution))
        push!(results, (; amplitude=target, result, solution, diagnostics))
        serialize(joinpath(directory,"ad_validation.jls"), results)
    end
    for candidate in (deuflhard, last(results))
        for (key, tolerance) in ((:energy, 1e-7), (:total_flux, 1e-5),
                                 (:omega_real, 1e-9), (:omega_imag, 1e-5))
            @test isapprox(getproperty(candidate.diagnostics,key),
                getproperty(direct.diagnostics,key); rtol=tolerance)
        end
        @test isapprox(candidate.solution.horizon_amplitudes, direct.solution.horizon_amplitudes; rtol=1e-6)
        @test isapprox(candidate.solution.outer_amplitudes, direct.solution.outer_amplitudes; rtol=1e-6)
    end
end
