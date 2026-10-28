"""
    TwoPhaseMatchedSolution(state, amplitudes, residual_norm, iterations)

Two-cloud solution of the radial matching equations at fixed prescribed
outer amplitudes, together with its convergence diagnostics.

- `state`: real vector packing the two complex frequencies, all scaled horizon
  coefficients, and the scaled free outer coefficients. All real parts precede
  all imaginary parts.
- `amplitudes`: the two prescribed outer coefficients, in physical units and
  cloud order. These are supplied separately from the unknowns in state.
- `residual_norm`: Euclidean norm of the real matching residual divided by
  coefficient_scale, using the same normalization as the Newton tolerance.
- `iterations`: number of Newton iterations used to reach the solution.

The state, amplitudes, and residual norm share a real numeric type T. Keep the
associated TwoPhaseProblem to interpret the coefficient ordering and
scale. Construction stores the supplied values; it does not solve the matching
equations or verify convergence.
"""
struct TwoPhaseMatchedSolution{T<:Real}
    state::Vector{T}
    amplitudes::Vector{T}
    residual_norm::T
    iterations::Int
end

"""
    two_phase_profile_admissible(packed_state, problem)

Ensure that a profile has sensible values for frequency and boundary coefficients.
"""
function two_phase_profile_admissible(packed_profile, problem::TwoPhaseProblem)
    all(isfinite, packed_profile) || return false
    frequencies = unpack_complex(packed_profile)[1:2]

    # Check that all (dimensionless, M*omega) frequencies are between 0 and the axion mass
    all(value -> 0 < real(value) < problem.alpha, frequencies) || return false
    component_real_frequencies = [real(sum(component.q .* frequencies))
        for component in problem.components]

    # Check that all component frequencies maintain their size relative to alpha
    return (component_real_frequencies .< problem.alpha) == problem.initial_bound
end

"""
    two_phase_real_matching_residual(packed_state, amplitudes, problem)

TODO: describe the real matching equations and residual scaling.
"""
function two_phase_real_matching_residual(packed_state, amplitudes, problem::TwoPhaseProblem)
    profile = unpack_two_phase_profile(packed_state, amplitudes, problem)
    residual = two_phase_matching_residual(profile, problem)
    return pack_complex(residual) ./ problem.coefficient_scale
end

"""
    correct_two_phase_profile(problem, amplitudes, packed_state_guess; kwargs...)

TODO: describe Newton correction, convergence tolerance, and returned packed state.
"""
function correct_two_phase_profile(problem::TwoPhaseProblem, amplitudes, packed_state_guess;
        tolerance=1e-9, max_iterations=12, progress_callback=nothing)
    residual = packed_state -> two_phase_real_matching_residual(packed_state, amplitudes, problem)
    admissible = packed_state -> two_phase_profile_admissible(packed_state, problem)
    solution = newton_solve(residual, packed_state_guess; tolerance, max_iterations,
        admissible, progress_callback, jacobian_backend=:forwarddiff)
    return (; state=solution.state, amplitudes=copy(amplitudes),
        residual_norm=norm(solution.residual), iterations=solution.iterations)
end

"""
    two_phase_trial_observables(packed_state, amplitudes, problem)

TODO: describe finite-domain charges, signed losses, and the real-carrier choice.
"""
function two_phase_trial_observables(packed_state, amplitudes, problem::TwoPhaseProblem)
    profile = unpack_two_phase_profile(packed_state, amplitudes, problem)
    radial_matching_result = two_phase_matching_residual(profile, problem; save_radial_result=true)
    # The matching frequencies are complex trial eigenvalues. The frozen-field
    # stress tensor uses their REAL parts, omitting slow profile time derivatives.
    # This explicit adapter does not establish the physical slow-growth closure.
    carrier_frequencies = real.(profile.frequencies)
    return two_phase_integrated_observables(radial_matching_result, problem.spin, problem.alpha,
        carrier_frequencies, problem.components, problem.cloud_m,
        problem.polar_values, problem.grid, problem.phase1_nodes, problem.phase2_nodes;
        potential=problem.potential)
end

"""
    two_phase_profile_response(problem, solution; tolerance=1e-9)

TODO: describe the E/J derivative rows and matching-response diagnostics.
"""
function two_phase_profile_response(problem::TwoPhaseProblem, solution; tolerance=1e-9)
    packed_state, amplitudes = solution.state, solution.amplitudes
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("require positive tolerance"))
    two_phase_profile_admissible(packed_state, problem) || error("inadmissible response state")
    residual = (trial_packed_state, trial_amplitudes) ->
        two_phase_real_matching_residual(trial_packed_state, trial_amplitudes, problem)
    residual_norm = norm(residual(packed_state, amplitudes))
    isfinite(residual_norm) && residual_norm <= tolerance ||
        error("profile must be converged before differentiation: residual=$residual_norm")
    charges = (trial_packed_state, trial_amplitudes) -> begin
        observables = two_phase_trial_observables(trial_packed_state, trial_amplitudes, problem)
        [observables.energy, observables.angular_momentum]
    end
    response = implicit_observable_derivatives(residual, charges, packed_state, amplitudes)
    observables = two_phase_trial_observables(packed_state, amplitudes, problem)
    return (; response..., observables, residual_norm)
end
