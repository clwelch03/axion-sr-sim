"""
    difference_jacobian(residual, state, increment)

Compute the production centered-difference Jacobian for comparison with AD.
"""
function difference_jacobian(residual, state, increment)
    matrix = zeros(length(residual(state)), length(state))
    for column in eachindex(state)
        displacement = increment * max(1, abs(state[column]))
        plus, minus = copy(state), copy(state)
        plus[column] += displacement
        minus[column] -= displacement
        matrix[:, column] = (residual(plus)-residual(minus))/(2displacement)
    end
    return matrix
end

"""
    correct_system(residual, jacobian, initial; tolerance=1e-9,
                   max_iterations=20, strategy=:residual, admissible, callback)

Correct a guess using residual backtracking or simplified-correction damping.

With `strategy=:deuflhard`, reuse the current factorization to measure the
simplified Newton correction at a trial. Require its norm to contract and
estimate damping from its departure from the linear model. This is a
Deuflhard-inspired monotonicity controller, not a full NLEQ-ERR implementation.
Always retain the same residual acceptance tolerance. Return failures as data
with the actual last iterate and iteration history for subsequent diagnosis.
"""
function correct_system(residual, jacobian, initial; tolerance=1e-9,
                        max_iterations=20, strategy=:residual,
                        admissible=state -> all(isfinite, state), callback=result -> nothing)
    state = copy(initial)
    history = NamedTuple[]
    residual_calls, jacobian_calls = 0, 0
    evaluate = function(state)
        residual_calls += 1
        residual(state)
    end
    mismatch = evaluate(state)
    status = "iteration_limit"
    for iteration in 0:max_iterations
        residual_norm = norm(mismatch)
        record = (; iteration, residual=residual_norm)
        push!(history, record)
        callback((; state=copy(state), record))
        println("iteration=", iteration, " residual=", residual_norm)
        flush(stdout)
        if !isfinite(residual_norm)
            status = "nonfinite_residual"
            break
        elseif residual_norm <= tolerance
            status = "converged"
            break
        elseif iteration == max_iterations
            break
        end
        matrix = jacobian(state)
        jacobian_calls += 1
        factorization = lu(matrix)
        correction = -(factorization \ mismatch)
        damping = 1.0
        accepted = false
        for attempt in 0:20
            trial = state + damping*correction
            if !admissible(trial)
                damping /= 2
                continue
            end
            trial_mismatch = evaluate(trial)
            if all(isfinite, trial_mismatch)
                if strategy === :residual
                    accepted = norm(trial_mismatch) < residual_norm
                elseif strategy === :deuflhard
                    simplified = -(factorization \ trial_mismatch)
                    accepted = norm(simplified) < norm(correction)
                    if !accepted
                        departure = norm(simplified-(1-damping)*correction)
                        estimate = departure == 0 ? damping/2 :
                            norm(correction)*damping^2/(2departure)
                        damping = min(damping/2, estimate)
                        continue
                    end
                else
                    error("Unknown correction strategy: $strategy")
                end
            end
            if accepted
                state, mismatch = trial, trial_mismatch
                println("accepted damping=", damping)
                break
            end
            damping /= 2
        end
        if !accepted
            status = "failed_descent"
            break
        end
    end
    return (; state, residual=norm(mismatch), status, history, residual_calls, jacobian_calls)
end
