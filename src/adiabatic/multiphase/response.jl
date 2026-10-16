"""
    two_phase_response_problem(reference_frequencies, reference_amplitudes,
        spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes, phase2_nodes; kwargs...)

TODO: describe fixed numerical settings, scale, and frequency-branch restrictions.
"""
function two_phase_response_problem(reference_frequencies, reference_amplitudes,
        spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes, phase2_nodes; potential=:cosine,
        coefficient_scale=minimum(reference_amplitudes), step=0.05,
        inner=-100.0, outer=100.0, match=5.0, horizon_order=3, outer_order=7)
    length(reference_frequencies) == length(reference_amplitudes) == length(cloud_m) == 2 ||
        throw(DimensionMismatch("require two independent clouds"))
    all(isfinite, (inner, match, outer, step)) && inner < match < outer && step > 0 ||
        throw(ArgumentError("require ordered endpoints and a positive finite radial step"))
    pairs = [(l, component.q) for component in components for l in component.l_values]
    fixed = [findfirst(pair -> pair[2] == q, pairs) for q in ((1, 0), (0, 1))]
    any(isnothing, fixed) && throw(ArgumentError("require both fundamental components"))
    free = setdiff(1:length(pairs), Int.(fixed))
    reference_state = complex.(vcat(reference_frequencies, zeros(2length(pairs)-2)))
    unpack_two_phase_trial(reference_state, reference_amplitudes, components; coefficient_scale)

    component_frequencies = [real(sum(component.q .* reference_frequencies))
        for component in components]
    initial_bound = component_frequencies .< alpha
    integration = (; step, inner, outer, match, horizon_order, outer_order)
    problem = (; spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes, phase2_nodes, potential, coefficient_scale, free,
        initial_bound, integration)
    two_phase_trial_admissible(pack_complex(reference_state), problem) ||
        throw(ArgumentError("reference frequencies are outside the supported boundary branches"))
    return problem
end

"""
    pack_two_phase_solution(solution, problem)

TODO: describe real packing and the fixed coefficient scale.
"""
function pack_two_phase_solution(solution, problem)
    scaled_coefficients = vcat(solution.horizon_coefficients,
        solution.outer_coefficients[problem.free]) ./ problem.coefficient_scale
    return pack_complex(vcat(solution.frequencies, scaled_coefficients))
end

"""
    two_phase_trial_admissible(state, problem)

TODO: describe the positive-frequency and bound/outgoing branch checks.
"""
function two_phase_trial_admissible(state, problem)
    all(isfinite, state) || return false
    frequencies = unpack_complex(state)[1:2]
    all(value -> 0 < real(value) < problem.alpha, frequencies) || return false
    component_frequencies = [real(sum(component.q .* frequencies))
        for component in problem.components]
    return all(value -> isfinite(value) && value > 0 && value != problem.alpha,
        component_frequencies) && (component_frequencies .< problem.alpha) == problem.initial_bound
end

"""
    two_phase_trial_matching(state, amplitudes, problem; save_profile=false)

TODO: describe reconstruction of physical trial coefficients and radial histories.
"""
function two_phase_trial_matching(state, amplitudes, problem; save_profile=false)
    trial = unpack_two_phase_trial(unpack_complex(state), amplitudes,
        problem.components; coefficient_scale=problem.coefficient_scale)
    return two_phase_matching_residual(trial..., problem.spin, problem.alpha,
        problem.components, problem.cloud_m, problem.couplings, problem.polar_values,
        problem.grid, problem.phase1_nodes, problem.phase2_nodes;
        potential=problem.potential, problem.integration..., save_profile)
end

"""
    two_phase_real_residual(state, amplitudes, problem)

TODO: describe the real matching equations and residual scaling.
"""
function two_phase_real_residual(state, amplitudes, problem)
    residual = two_phase_trial_matching(state, amplitudes, problem)
    return pack_complex(residual) ./ problem.coefficient_scale
end

"""
    correct_two_phase_profile(problem, amplitudes, state_guess; kwargs...)

TODO: describe Newton correction, convergence tolerance, and returned state.
"""
function correct_two_phase_profile(problem, amplitudes, state_guess;
        tolerance=1e-9, max_iterations=12, progress_callback=nothing)
    residual = state -> two_phase_real_residual(state, amplitudes, problem)
    admissible = state -> two_phase_trial_admissible(state, problem)
    solution = newton_solve(residual, state_guess; tolerance, max_iterations,
        admissible, progress_callback, jacobian_backend=:forwarddiff)
    return (; state=solution.state, amplitudes=copy(amplitudes),
        residual_norm=norm(solution.residual), iterations=solution.iterations)
end

"""
    two_phase_trial_observables(state, amplitudes, problem)

TODO: describe finite-domain charges, signed losses, and the real-carrier choice.
"""
function two_phase_trial_observables(state, amplitudes, problem)
    profiles = two_phase_trial_matching(state, amplitudes, problem; save_profile=true)
    # The matching frequencies are complex trial eigenvalues. The frozen-field
    # stress tensor uses their REAL parts, omitting slow profile time derivatives.
    # This explicit adapter does not establish the physical slow-growth closure.
    carrier_frequencies = real.(unpack_complex(state)[1:2])
    return two_phase_integrated_observables(profiles, problem.spin, problem.alpha,
        carrier_frequencies, problem.components, problem.cloud_m,
        problem.polar_values, problem.grid, problem.phase1_nodes, problem.phase2_nodes;
        potential=problem.potential)
end

"""
    two_phase_profile_response(problem, solution; tolerance=1e-9)

TODO: describe the E/J derivative rows and matching-response diagnostics.
"""
function two_phase_profile_response(problem, solution; tolerance=1e-9)
    state, amplitudes = solution.state, solution.amplitudes
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("require positive tolerance"))
    two_phase_trial_admissible(state, problem) || error("inadmissible response state")
    residual = (trial_state, trial_amplitudes) ->
        two_phase_real_residual(trial_state, trial_amplitudes, problem)
    residual_norm = norm(residual(state, amplitudes))
    isfinite(residual_norm) && residual_norm <= tolerance ||
        error("profile must be converged before differentiation: residual=$residual_norm")
    charges = (trial_state, trial_amplitudes) -> begin
        observables = two_phase_trial_observables(trial_state, trial_amplitudes, problem)
        [observables.energy, observables.angular_momentum]
    end
    response = implicit_observable_derivatives(residual, charges, state, amplitudes)
    observables = two_phase_trial_observables(state, amplitudes, problem)
    return (; response..., observables, residual_norm)
end
