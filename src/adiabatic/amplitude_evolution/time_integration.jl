"""
    rk4_amplitude_step(rate, amplitude, time_step)

Advance the autonomous equation dA/dt = rate(A) by one RK4 step.
Time is measured in t/M when using the tabulated amplitude_rate.

This function does not choose a safe step size or enforce branch boundaries.
"""
function rk4_amplitude_step(rate, amplitude, time_step)
    isfinite(amplitude) ||
        throw(ArgumentError("require finite amplitude"))
    isfinite(time_step) && time_step > 0 ||
        throw(ArgumentError("require a positive finite time step"))

    first_rate = rate(amplitude)
    second_rate = rate(amplitude + time_step * first_rate / 2)
    third_rate = rate(amplitude + time_step * second_rate / 2)
    fourth_rate = rate(amplitude + time_step * third_rate)

    next_amplitude = amplitude + time_step *
        (first_rate + 2second_rate + 2third_rate + fourth_rate) / 6

    all(isfinite, (first_rate, second_rate, third_rate,
                   fourth_rate, next_amplitude)) ||
        error("nonfinite amplitude or rate during RK4 step")

    return next_amplitude
end

function checked_amplitude_rate(amplitude, branch)
    interval = branch_interval(amplitude, branch.amplitudes)
    left = branch.samples[interval.left_index]
    right = branch.samples[interval.right_index]
    derivative = interpolate_branch_energy(
        amplitude, left, right).energy_derivative
    isfinite(derivative) && derivative > 0 ||
        throw(DomainError(
            amplitude, "require positive finite energy derivative"))
    return amplitude_rate_on_interval(amplitude, left, right)
end

"""
    evolve_amplitude(branch, initial_amplitude; final_time, time_step)

Integrate dA0/d(t/M) on a branch with positive energy derivative.
Use uniform RK4 steps no larger than time_step, ending at final_time.
Reject evaluations outside the table or at nonpositive energy derivative.
Check timestep convergence separately; time_step is not an error tolerance.
"""
function evolve_amplitude(branch, initial_amplitude; final_time, time_step)
    isfinite(final_time) && final_time >= 0 ||
        throw(ArgumentError("require finite nonnegative final time"))
    isfinite(time_step) && time_step > 0 ||
        throw(ArgumentError("require positive finite time step"))
    length(branch.amplitudes) == length(branch.samples) ||
        throw(DimensionMismatch("require one sample per amplitude"))
    all(sample.amplitude == amplitude for
        (sample, amplitude) in zip(branch.samples, branch.amplitudes)) ||
        throw(ArgumentError("sample amplitudes do not match the grid"))

    rate = amplitude -> checked_amplitude_rate(amplitude, branch)
    isfinite(rate(initial_amplitude)) || error("nonfinite initial rate")
    if iszero(final_time)
        return (; times=[0.0], amplitudes=[Float64(initial_amplitude)])
    end

    step_count = ceil(Int, final_time / time_step)
    times = collect(range(0.0, Float64(final_time); length=step_count + 1))
    amplitudes = Vector{Float64}(undef, step_count + 1)
    amplitudes[1] = initial_amplitude

    for index in 1:step_count
        next_amplitude = rk4_amplitude_step(
            rate, amplitudes[index], times[index+1] - times[index])
        isfinite(rate(next_amplitude)) || error("nonfinite accepted rate")
        amplitudes[index+1] = next_amplitude
    end

    return (; times, amplitudes)
end
