using ForwardDiff

"""
    implicit_branch_derivative(residual, state, amplitude)

Compute the state derivative along a branch satisfying residual(state, amplitude) = 0.

Supply a real state vector at a converged solution and a residual function
returning a real vector of the same length. The residual must support
ForwardDiff differentiation in both arguments.

Evaluate the state Jacobian and the partial amplitude derivative, then
solve the differentiated matching equations for `dstate/damplitude`.
Require a nonsingular state Jacobian.

Return the real derivative vector in the same ordering as `state`.
The caller is responsible for checking convergence of the supplied solution.
"""
function implicit_branch_derivative(
    residual, state::AbstractVector{<:Real}, amplitude::Real,
)
    state_jacobian = ForwardDiff.jacobian(
        trial_state -> residual(trial_state, amplitude),
        state,
    )

    amplitude_derivative = ForwardDiff.jacobian(
        trial_amplitude -> residual(state, trial_amplitude[1]),
        [amplitude],
    )[:, 1]

    return -(state_jacobian \ amplitude_derivative)
end

"""
    implicit_energy_derivative(residual, energy, state, amplitude)

Compute the total energy derivative along a stationary solution branch.

Supply the real residual and converged state required by
[`implicit_branch_derivative`](@ref), together with a function
`energy(state, amplitude)` returning the real scalar energy. Both functions
must support ForwardDiff differentiation. Keep the physical parameters and
spatial discretization fixed while differentiating.

Solve for the state derivative, then differentiate energy along that direction,
including its explicit amplitude dependence. Return a named tuple containing
`state_derivative` and `energy_derivative`. The latter preserves the supplied
energy normalization; input `E/Fa²` gives `d(E/Fa²)/dA0`.

This differentiates the discretized branch and does not estimate its spatial
discretization error. The caller must supply a converged solution with a
nonsingular state Jacobian.

See also [`implicit_branch_derivative`](@ref).
"""
function implicit_energy_derivative(
    residual, energy, state::AbstractVector{<:Real}, amplitude::Real,
)
    result = implicit_observable_derivatives(
        residual,
        (trial_state, trial_amplitude) -> [energy(trial_state, trial_amplitude)],
        state,
        amplitude,
    )

    return (;
        state_derivative=result.state_derivative,
        energy_derivative=only(result.observable_derivatives),
    )
end

"""
    implicit_observable_derivatives(residual, observables, state, amplitude)

Compute total amplitude derivatives along a stationary solution branch.

`observables(state, amplitude)` must return a real vector and support
ForwardDiff, for example `[energy, horizon_flux, outer_flux]`.

Require a converged state and a nonsingular residual Jacobian. Keep
physical parameters and spatial settings fixed while differentiating.

Return the state derivative and observable derivatives in the supplied
order. This does not estimate spatial discretization error.
"""
function implicit_observable_derivatives(
    residual, observables, state::AbstractVector{<:Real}, amplitude::Real,
)
    state_derivative = implicit_branch_derivative(residual, state, amplitude)

    observable_derivatives = ForwardDiff.jacobian(
        offsets -> observables(
            state .+ offsets[1] .* state_derivative,
            amplitude + offsets[1],
        ),
        [zero(amplitude)],
    )[:, 1]

    return (; state_derivative, observable_derivatives)
end