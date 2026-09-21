using LinearAlgebra: norm


"""
    newton_solve(residual, initial; <keyword arguments>)

Solve `residual(state) = 0` using damped Newton steps and centered differences.

Leave `initial` unchanged and return a named tuple with `state`, `residual`,
and `iterations`. Real initial values produce a real Jacobian; complex initial
values require a residual analytic in its complex arguments. Throw an error on
nonfinite residuals, inadmissible differences, failed descent, or iteration exhaustion.

# Arguments
- `tolerance=1e-10`: positive finite upper bound on the residual norm.
- `max_iterations=15`: nonnegative integer maximum number of Newton steps.
- `difference_step=1e-6`: positive finite difference scale, multiplied by
  `max(1, abs(state[j]))` for each component.
- `admissible=state -> true`: predicate restricting guesses and trial states;
  both points of each centered difference must also satisfy it.

See also [`real_newton`](@ref), [`complex_newton`](@ref).
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


"""
    pack_complex(values)

Pack complex values into a vector of real parts followed by imaginary parts.

See also [`unpack_complex`](@ref).

# Examples
```jldoctest
julia> pack_complex([1 + 2im, 3 + 4im]) == [1, 3, 2, 4]
true
```
"""
pack_complex(values) = vcat(real.(values), imag.(values))

"""
    unpack_complex(values::AbstractVector{<:Real})

Recover a complex vector from equally sized real and imaginary halves.

Require even length, and interpret `values` in the layout of [`pack_complex`](@ref).
"""
function unpack_complex(values::AbstractVector{<:Real})
    iseven(length(values)) ||
        throw(DimensionMismatch("expected equally sized real and imaginary halves"))

    count = length(values) ÷ 2
    return complex.(values[1:count], values[count+1:end])
end

"""
    real_newton(residual, initial; admissible=state -> true, kwargs...)

Solve a possibly nonanalytic complex system using a real Jacobian.

Pass complex vectors to `residual` and `admissible`, while packing real and
imaginary parts separately for [`newton_solve`](@ref). Forward `kwargs...`
(`tolerance`, `max_iterations`, and `difference_step`) to that solver with its
defaults. Return a named tuple with complex `state`, complex `residual`, and
`iterations`; leave `initial` unchanged.

See also [`complex_newton`](@ref).
"""
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


"""
    complex_newton(residual, initial; kwargs...)

Solve an analytic complex system using a complex Jacobian.

Convert `initial` to complex values and forward keyword arguments to
[`newton_solve`](@ref), including its defaults for `tolerance`, `max_iterations`,
`difference_step`, and `admissible`. Return its named tuple containing `state`,
`residual`, and `iterations`; leave `initial` unchanged.

See also [`real_newton`](@ref) for residuals involving complex conjugates.
"""
function complex_newton(residual, initial; kwargs...)
    return newton_solve(residual, complex.(initial); kwargs...)
end


"""
    solve_linear_mode(frequency_guess, spin, alpha, l_values, m, couplings;
                      <keyword arguments>)

Find a linear bound mode with the first outer amplitude fixed to one.

Put the normalization mode first in `l_values`; use `[1, 3, 5]` and `m = 1`
for the fundamental cloud. Solve for the frequency, all horizon amplitudes,
and the remaining outer amplitudes in units `M = 1`. Require the real frequency
to lie strictly between zero and `alpha` during Newton iteration.

Return a named tuple with `frequency`, `horizon_amplitudes`, `outer_amplitudes`,
`residual`, `residual_norm`, and `iterations`. A small mismatch at the selected
settings does not establish endpoint or step convergence, especially for the
small imaginary frequency.

# Arguments
- `horizon_guess=ones(length(l_values))`: initial horizon amplitudes in `l_values` order.
- `outer_tail_guess=zeros(length(l_values)-1)`: initial free outer amplitudes.
- `tolerance=1e-10`: positive finite matching-residual norm tolerance.
- `max_iterations=15`: nonnegative integer Newton iteration limit.
- `difference_step=1e-6`: positive finite relative scale for centered differences.
- `integration_options...`: keywords for [`matching_residual`](@ref): `step=0.05`,
  `inner=-100.0`, `outer=100.0`, and `match=5.0`.

See also [`solve_nonlinear_mode`](@ref).
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

"""
    solve_nonlinear_mode(frequency_guess, amplitude, spin, alpha, modes, couplings,
                         polar_values, azimuthal_values, grid; <keyword arguments>)

Solve nonlinear matching at the fixed positive outer amplitude `amplitude = A0`.

Put the fundamental normalization mode first in `modes`. Supply physical
boundary amplitudes as guesses; Newton scales them by `A0` and applies its
tolerance to `residual/A0`. Require positive finite real `amplitude` and keep
the fundamental frequency on the bound branch in units `M = 1`.

Return a named tuple with `frequency`, physical `horizon_amplitudes` and
`outer_amplitudes`, physical `residual` and `residual_norm`, `scaled_residual_norm`,
and `iterations`. The first outer amplitude remains `A0`. Check step, endpoint,
series-order, and angular convergence separately from the matching tolerance.

# Arguments
- `horizon_guess`: required vector of one physical amplitude per mode.
- `outer_tail_guess`: required vector of physical outer amplitudes for `modes[2:end]`.
- `fundamental_m=1`: positive integer azimuthal number of the fundamental.
- `tolerance=1e-10`: positive finite tolerance on the scaled residual norm.
- `max_iterations=15`: nonnegative integer Newton iteration limit.
- `difference_step=1e-6`: positive finite relative scale for centered differences.
- `integration_options...`: keywords for [`nonlinear_matching_residual`](@ref):
  `step=0.05`, `inner=-100.0`, `outer=100.0`, `match=5.0`, `horizon_order=3`,
  and `outer_order=7`.

See also [`solve_linear_mode`](@ref), [`real_newton`](@ref).
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
