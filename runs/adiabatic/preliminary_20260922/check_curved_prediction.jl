# Diagnostic only: compare a curved predictor with the previous secant trial.
include(joinpath(@__DIR__, "check_residual_floor.jl"))

function curved_prediction(continuation, target)
    last_amplitude = last(continuation).amplitude
    # Avoid estimating curvature from the tightly clustered final points.
    indices = [findlast(entry -> entry.amplitude <= last_amplitude-offset, continuation)
        for offset in (0.15, 0.10, 0.05, 0.0)]
    entries = continuation[indices]
    nodes = [entry.amplitude for entry in entries]
    length(unique(nodes)) == 4 || error("require four distinct predictor nodes")
    weights = [prod((target-nodes[other])/(nodes[index]-nodes[other])
        for other in eachindex(nodes) if other != index) for index in eachindex(nodes)]
    frequency = sum(weight*entry.solution.frequency for (weight, entry) in zip(weights, entries))
    horizon_amplitudes = sum(weight .* entry.solution.horizon_amplitudes
        for (weight, entry) in zip(weights, entries))
    outer_amplitudes = sum(weight .* entry.solution.outer_amplitudes
        for (weight, entry) in zip(weights, entries))
    outer_amplitudes[1] = target
    println("CUBIC target=", target, " nodes=", nodes, " weights=", weights)
    flush(stdout)
    return (; frequency, horizon_amplitudes, outer_amplitudes)
end

target = 10.0
trial = solve_trial(saved.continuation, target, 5.0, 1e-10;
    prediction=curved_prediction(saved.continuation, target))
serialize(joinpath(@__DIR__, "curved_prediction_diagnosis.jls"), (; seed_path, trial))
