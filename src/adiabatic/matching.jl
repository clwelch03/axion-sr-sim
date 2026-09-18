using LinearAlgebra: norm


"""Solve F(state)=0 using damped Newton steps and centered differences.

Real initial values produce a real Jacobian, allowing general real systems.
Complex initial values require F to be analytic in its complex arguments.
"""
function newton_solve(residual, initial; tolerance=1e-10, max_iterations=15,
                        difference_step=1e-6, admissible=state -> true)
    isfinite(tolerance) && tolerance > 0 && isfinite(difference_step) && difference_step > 0 &&
        max_iterations isa Integer && max_iterations >= 0 ||
        throw(ArgumentError("require positive tolerances and nonnegative iteration limit"))
    state = float.(initial)
    admissible(state) || throw(ArgumentError("initial guess is outside the allowed domain"))
    mismatch = residual(state)
    for iteration in 0:max_iterations
        mismatch_norm = norm(mismatch)
        isfinite(mismatch_norm) || error("nonfinite Newton residual")
        mismatch_norm <= tolerance && return (; state, residual=mismatch, iterations=iteration)
        iteration == max_iterations && error("Newton iteration limit reached; residual=$mismatch_norm")

        # J[i,j] measures how matching condition i changes with unknown j.
        jacobian = Matrix{eltype(state)}(undef, length(mismatch), length(state))
        for column in eachindex(state)
            displacement = difference_step * max(1, abs(state[column]))
            plus, minus = copy(state), copy(state)
            plus[column] += displacement
            minus[column] -= displacement
            admissible(plus) && admissible(minus) ||
                error("Jacobian difference crosses the domain boundary; reduce difference_step")
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2displacement)
        end
        # F(z+dz) ~ F(z)+J*dz: solve for the correction, without inverting J.
        correction = -(jacobian \ mismatch)
        accepted = false
        for reduction in 0:20
            trial = state + 2.0^(-reduction) * correction
            admissible(trial) || continue
            trial_mismatch = residual(trial)
            if norm(trial_mismatch) < mismatch_norm
                state, mismatch = trial, trial_mismatch
                accepted = true
                break
            end
        end
        accepted || error("Newton step could not reduce the residual ($mismatch_norm)")
    end
end


"""Store all real parts first, followed by all imaginary parts."""
pack_complex(values) = vcat(real.(values), imag.(values))

"""Recover complex values from [real parts; imaginary parts]."""
function unpack_complex(values::AbstractVector{<:Real})
    iseven(length(values)) ||
        throw(DimensionMismatch("expected equally sized real and imaginary halves"))

    count = length(values) ÷ 2
    return complex.(values[1:count], values[count+1:end])
end

"""Solve a possibly nonanalytic complex system using a real Jacobian."""
function real_newton(residual, initial;
                     admissible=state -> true, kwargs...)
    # Physics functions receive complex values; Newton receives real vectors.
    real_residual = packed -> pack_complex(
        residual(unpack_complex(packed))
    )
    real_admissible = packed -> admissible(unpack_complex(packed))

    solution = newton_solve(
        real_residual, pack_complex(initial);
        admissible=real_admissible, kwargs...,
    )

    # Return the same complex-valued interface as complex_newton.
    return (
        state=unpack_complex(solution.state),
        residual=unpack_complex(solution.residual),
        iterations=solution.iterations,
    )
end


"""Solve an analytic complex system using a complex Jacobian."""
function complex_newton(residual, initial; kwargs...)
    return newton_solve(residual, complex.(initial); kwargs...)
end


"""Find a linear bound mode, fixing outer_amplitudes[1]=1.

Use l_values=[1,3,5], m=1 for the fundamental cloud. The fixed coefficient
belongs to the FIRST supplied l, so put the mode used for normalization first.
For N modes the 2N complex unknowns are [frequency; N horizon amplitudes;
N-1 outer amplitudes]. Matching supplies 2N complex equations for them.

step, inner, outer, and match are forwarded to matching_residual. A small
residual establishes matching at these numerical settings; endpoint and step
convergence are separate checks, especially for the small imaginary frequency.
"""
function solve_linear_mode(frequency_guess, spin, alpha, l_values, m, couplings;
                           horizon_guess=ones(length(l_values)),
                           outer_tail_guess=zeros(length(l_values) - 1),
                           tolerance=1e-10, max_iterations=15,
                           difference_step=1e-6, integration_options...)
    count = length(l_values)
    count > 0 && length(horizon_guess) == count && length(outer_tail_guess) == count - 1 ||
        throw(DimensionMismatch("require N horizon and N-1 free outer amplitudes"))

    initial = vcat(complex(frequency_guess), complex.(horizon_guess), complex.(outer_tail_guess))
    residual = state -> matching_residual(
        state[1], state[2:count+1], vcat(one(state[1]), state[count+2:end]),
        spin, alpha, l_values, m, couplings; integration_options...)

    # Keep Newton trials on the positive-frequency bound branch.
    admissible = state -> all(isfinite, state) && 0 < real(state[1]) < alpha
    solution = complex_newton(residual, initial;
        tolerance, max_iterations, difference_step, admissible)

    return (; frequency=solution.state[1],
              horizon_amplitudes=solution.state[2:count+1],
              outer_amplitudes=vcat(one(solution.state[1]), solution.state[count+2:end]),
              residual=solution.residual, residual_norm=norm(solution.residual),
              iterations=solution.iterations)
end

"""Solve the nonlinear matching problem at fixed positive outer amplitude A0.

Put the fundamental normalization mode first in modes.
horizon_guess and outer_tail_guess contain physical boundary amplitudes.
Newton uses amplitudes scaled by A0; tolerance applies to residual/A0.
"""
function solve_nonlinear_mode(
    frequency_guess, amplitude, spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    horizon_guess, outer_tail_guess,
    fundamental_m=1, tolerance=1e-10,
    max_iterations=15, difference_step=1e-6,
    integration_options...,
)
    amplitude isa Real && isfinite(amplitude) && amplitude > 0 ||
        throw(ArgumentError("require positive finite A0"))

    count = length(modes)
    count > 0 && length(horizon_guess) == count &&
        length(outer_tail_guess) == count - 1 ||
        throw(DimensionMismatch("require N horizon and N-1 outer guesses"))
    modes[1][2] == fundamental_m ||
        throw(ArgumentError("put the fundamental normalization mode first"))

    # Unknowns: frequency, N scaled horizon amplitudes, N-1 scaled outer ones.
    # The first scaled outer amplitude is fixed to 1, so its physical value is A0.
    initial = vcat(
        complex(frequency_guess),
        complex.(horizon_guess ./ amplitude),
        complex.(outer_tail_guess ./ amplitude),
    )

    residual = state -> nonlinear_matching_residual(
        state[1],
        amplitude .* state[2:count+1],
        amplitude .* vcat(one(state[1]), state[count+2:end]),
        spin, alpha, modes, couplings,
        polar_values, azimuthal_values, grid;
        fundamental_m, integration_options...,
    ) ./ amplitude

    # The fundamental remains bound. Higher harmonics use their own
    # boundary branches, as selected by boundary_state.
    admissible = state ->
        all(isfinite, state) && 0 < real(state[1]) < alpha

    # Reconstruction involves conjugates, so use the real Jacobian.
    solution = real_newton(
        residual, initial;
        tolerance, max_iterations, difference_step, admissible,
    )

    physical_residual = amplitude .* solution.residual
    return (
        frequency=solution.state[1],
        horizon_amplitudes=amplitude .* solution.state[2:count+1],
        outer_amplitudes=amplitude .* vcat(
            one(solution.state[1]), solution.state[count+2:end],
        ),
        residual=physical_residual,
        residual_norm=norm(physical_residual),
        scaled_residual_norm=norm(solution.residual),
        iterations=solution.iterations,
    )
end