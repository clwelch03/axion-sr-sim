# Test whether the locally improved Jacobian permits useful amplitude steps.
include(joinpath(@__DIR__, "check_residual_floor.jl"))

function check_finer_correction()
    continuation = copy(saved.continuation)
    trials = NamedTuple[]
    for target in (10.0, 10.5)
        trial = solve_trial(continuation, target, 5.0, 1e-10)
        push!(trials, trial)
        serialize(joinpath(@__DIR__, "finer_correction_diagnosis.jls"), (; seed_path, trials))
        isnothing(trial.solution) && break
        push!(continuation, (; amplitude=target, solution=trial.solution))
    end
end

check_finer_correction()
