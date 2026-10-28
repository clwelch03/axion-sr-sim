using LinearAlgebra: norm
using ForwardDiff

"""
    residual_jacobian(residual, state; backend=:finite,
                      difference_step=1e-6, admissible=state -> true)

Compute the Jacobian of `residual` at `state`.

Select centered differences with `:finite` or automatic differentiation with
`:forwarddiff`. ForwardDiff requires a real state vector.

See also [`centered_difference_jacobian`](@ref), [`newton_solve`](@ref).
"""
function residual_jacobian(residual, state;
        backend=:finite, chunk_size=nothing,
        difference_step=1e-6, admissible=state -> true)

    if chunk_size !== nothing && !(backend === :forwarddiff &&
            chunk_size isa Int && 1 <= chunk_size <= length(state))
        throw(ArgumentError(
            "chunk_size requires :forwarddiff and an integer from 1 to length(state)"))
    end

    if backend === :forwarddiff
        config = if chunk_size === nothing
            ForwardDiff.JacobianConfig(residual, state)
        else
            ForwardDiff.JacobianConfig(
                residual, state, ForwardDiff.Chunk{chunk_size}())
        end
        return ForwardDiff.jacobian(residual, state, config)
    elseif backend === :finite
        return centered_difference_jacobian(
            residual, state; difference_step, admissible)
    end

    throw(ArgumentError("unknown Jacobian backend: $backend"))
end


"""
    centered_difference_jacobian(residual, state; difference_step=1e-6,
                                 admissible=state -> true)

Approximate each Jacobian column by centered differences.

Perturb component `j` by `difference_step * max(1, abs(state[j]))`.
Require both perturbed states to satisfy `admissible`. Leave `state` unchanged.

See also [`residual_jacobian`](@ref).
"""
function centered_difference_jacobian(
    residual,
    state;
    difference_step=1e-6,
    admissible=state -> true,
)
    columns = map(eachindex(state)) do index
        displacement = difference_step * max(1, abs(state[index]))

        forward_state = copy(state)
        backward_state = copy(state)
        forward_state[index] += displacement
        backward_state[index] -= displacement

        if !admissible(forward_state) || !admissible(backward_state)
            error(
                "Jacobian difference crosses the domain boundary; " *
                "reduce difference_step"
            )
        end

        forward_residual = residual(forward_state)
        backward_residual = residual(backward_state)

        return (forward_residual - backward_residual) / (2 * displacement)
    end

    return hcat(columns...)
end

"""
    newton_solve(residual, initial; <keyword arguments>)

Solve `residual(state) = 0` using damped Newton steps and centered differences.

Leave `initial` unchanged and return a named tuple with `state`, `residual`,
and `iterations`. Real initial values produce a real Jacobian; complex initial
values require a residual analytic in its complex arguments. Throw an error on
nonfinite residuals, inadmissible differences, failed descent, or iteration exhaustion.
With `max_iterations=nothing`, continue until the residual meets `tolerance`,
unless a correction fails or the residual stagnates.

# Arguments
- `tolerance=1e-10`: positive finite upper bound on the residual norm.
- `max_iterations=15`: nonnegative integer step limit, or `nothing` for no fixed limit.
- `difference_step=1e-6`: positive finite difference scale, multiplied by
  `max(1, abs(state[j]))` for each component.
- `admissible=state -> true`: predicate restricting guesses and trial states;
  both points of each centered difference must also satisfy it.
- `progress_callback=nothing`: optional function called with `(iteration, residual_norm)`
  at the initial state and after every accepted correction.
- `stagnation_limit=8`: positive number of consecutive insufficient improvements
  allowed before reporting stagnation.
- `minimum_progress=1e-3`: minimum fractional residual decrease per correction
  that resets the stagnation counter; require a value strictly between zero and one.

See also [`real_newton`](@ref), [`complex_newton`](@ref).
"""
function newton_solve(
    residual,
    initial;
    tolerance=1e-10,
    max_iterations=15,
    jacobian_backend=:finite,
    difference_step=1e-6,
    admissible=state -> true,
    progress_callback=nothing,
    stagnation_limit::Int=8,
    minimum_progress=1e-3,
    jacobian_chunk_size=nothing,
    accepted_state_callback=nothing,
)
    # Validate solver settings before evaluating the residual.
    if !isfinite(tolerance) || tolerance <= 0
        throw(ArgumentError("tolerance must be positive and finite"))
    end

    if !isnothing(max_iterations)
        if !(max_iterations isa Integer) || max_iterations < 0
            throw(ArgumentError(
                "max_iterations must be a nonnegative integer or nothing"
            ))
        end
    end
    if !isfinite(minimum_progress) || !(0 < minimum_progress < 1)
        throw(ArgumentError("minimum_progress must lie strictly between 0 and 1"))
    end

    if !(jacobian_backend in (:finite, :forwarddiff))
        throw(ArgumentError(
            "jacobian_backend must be :finite or :forwarddiff"
        ))
    end

    stagnation_limit > 0 || throw(ArgumentError("stagnation_limit must be positive"))

    # Work on a floating-point copy, leaving the caller's initial guess unchanged.
    state = float.(initial)

    if jacobian_backend === :forwarddiff && !(eltype(state) <: Real)
        throw(ArgumentError(
            "ForwardDiff requires the real-packed Newton system; use real_newton"
        ))
    end

    if !admissible(state)
        throw(ArgumentError("initial guess is outside the allowed domain"))
    end

    mismatch = residual(state)
    iteration = 0
    stagnant_steps = 0

    while true
        # 1. Check the current iterate before attempting another correction.
        mismatch_norm = norm(mismatch)

        if !isfinite(mismatch_norm)
            error("nonfinite Newton residual")
        end
        
        if !isnothing(progress_callback)
            progress_callback(iteration, mismatch_norm)
        end

        if !isnothing(accepted_state_callback)
            accepted_state_callback(iteration, copy(state), copy(mismatch))
        end

        if mismatch_norm <= tolerance
            return (
                state=state,
                residual=mismatch,
                iterations=iteration,
            )
        end

        if iteration == max_iterations
            error("Newton iteration limit reached; residual=$mismatch_norm")
        end

        if stagnant_steps >= stagnation_limit
            error(
                "Newton residual stagnated after $iteration corrections; " *
                "residual=$mismatch_norm"
            )
        end

        # 2. Linearize the residual around the current state.
        # J[i, j] measures how residual component i changes with unknown j.
        jacobian = residual_jacobian(
            residual, state;
            backend=jacobian_backend,
            chunk_size=jacobian_chunk_size,
            difference_step,
            admissible,
        )

        # 3. Solve J * correction = -mismatch.
        # Backslash solves the linear system without explicitly inverting J.
        correction = -(jacobian \ mismatch)

        # 4. Try the full Newton correction first.
        # If it leaves the allowed domain or fails to improve the residual,
        # repeatedly halve it. Only an accepted trial replaces the current state.
        accepted = false

        for reduction in 0:20
            damping = 2.0^(-reduction)
            trial_state = state + damping * correction

            if !admissible(trial_state)
                continue
            end

            trial_mismatch = residual(trial_state)

            if norm(trial_mismatch) < mismatch_norm
                state = trial_state
                mismatch = trial_mismatch
                accepted = true
                break
            end
        end

        if !accepted
            error("Newton step could not reduce the residual ($mismatch_norm)")
        end

        # 5. Track corrections that improve the residual only marginally.
        # A sufficiently large improvement resets the stagnation counter.
        required_residual = (1 - minimum_progress) * mismatch_norm

        if norm(mismatch) > required_residual
            stagnant_steps += 1
        else
            stagnant_steps = 0
        end

        iteration += 1
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
(`tolerance`, `max_iterations`, `difference_step`, and progress/stagnation
options) to that solver with its defaults. Return a named tuple with complex `state`, complex `residual`, and
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
