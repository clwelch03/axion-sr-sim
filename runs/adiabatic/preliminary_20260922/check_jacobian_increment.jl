# Separate local derivative study at the same preserved high-amplitude seed.
include(joinpath(@__DIR__, "check_residual_floor.jl"))

function check_jacobian_increment()
    previous = last(saved.continuation)
    amplitude, solution = previous.amplitude, previous.solution
    state = pack_complex(vcat(solution.frequency,
        solution.horizon_amplitudes ./ amplitude,
        solution.outer_amplitudes[2:end] ./ amplitude))
    residual = state -> packed_residual(state, amplitude, options)
    baseline = residual(state)
    results = NamedTuple[]
    for increment in (1e-8, 1e-9, 1e-10, 1e-11)
        jacobian = Matrix{Float64}(undef, length(baseline), length(state))
        for column in eachindex(state)
            displacement = increment * max(1, abs(state[column]))
            plus, minus = copy(state), copy(state)
            plus[column] += displacement
            minus[column] -= displacement
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2displacement)
        end
        correction = -(jacobian \ baseline)
        after_correction = residual(state + correction)
        column_norms = [norm(jacobian[:, column]) for column in eachindex(state)]
        result = (; increment, jacobian, correction, after_correction,
            condition=cond(jacobian), equilibrated_condition=cond(jacobian ./ transpose(column_norms)))
        push!(results, result)
        println("JACOBIAN increment=", increment, " condition=", result.condition,
            " equilibrated_condition=", result.equilibrated_condition,
            " initial_residual=", norm(baseline), " corrected_residual=", norm(after_correction))
        flush(stdout)
        serialize(joinpath(@__DIR__, "jacobian_increment_diagnosis.jls"), (; state, baseline, results))
    end
    reference = last(results).jacobian
    for result in results[1:end-1]
        column_errors = [norm(result.jacobian[:, column] - reference[:, column]) /
            norm(reference[:, column]) for column in eachindex(state)]
        println("COLUMN_ERRORS increment=", result.increment, " relative=", column_errors)
    end
end

check_jacobian_increment()
