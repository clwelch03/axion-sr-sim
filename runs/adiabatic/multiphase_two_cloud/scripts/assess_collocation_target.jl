"""
    assess_two_phase_collocation(result, problem::TwoPhaseProblem;
        sample_step=0.25, parameter_tolerance=1e-10)

Check a collocation result and return physical radial samples and diagnostics.

`result` is returned by `solve_two_phase_collocation`; supply its associated
`problem`. Evaluate both endpoint residuals, parameter drift across mesh nodes,
and the scaled differential defect at two independent points per interval.
Sample physical complex `[R; Q]` on an increasing grid covering both endpoints,
with spacing no greater than the positive finite `sample_step`.

Return `(; passed, profile, positions, states, boundary_residual,
parameter_drift, defect_by_state, absolute_defect_by_state, admissible)`.
`passed` requires library success, finite states and diagnostics, the intended
frequency branches, boundary residual within `result.settings.nonlinear_tolerance`,
maximum field defect within `result.settings.tolerance`, and parameter drift
within the positive finite `parameter_tolerance`. Check spatial refinement
separately before interpreting curve differences. Input arrays are read only.
"""
function assess_two_phase_collocation(result, problem::TwoPhaseProblem;
        sample_step=0.25, parameter_tolerance=1e-10)
    all(value -> isfinite(value) && value > 0, (sample_step, parameter_tolerance)) ||
        throw(ArgumentError("require positive finite sampling and parameter tolerances"))
    (; solution, amplitudes, state_scales, packed_reference, parameter_scales) = result
    field_count = length(state_scales)
    parameter_indices = (field_count + 1):length(first(solution.u))
    reference_offsets = view(first(solution.u), parameter_indices)
    parameter_drift = maximum(norm(view(state, parameter_indices) - reference_offsets, Inf)
        for state in solution.u)
    callbacks = two_phase_collocation_callbacks(problem, amplitudes;
        state_scales, packed_reference, parameter_scales)

    # Endpoint constraints use the current parameters at their respective ends.
    left_residual, right_residual = zeros(field_count), zeros(field_count)
    callbacks.horizon!(left_residual, first(solution.u), nothing)
    callbacks.outer!(right_residual, last(solution.u), nothing)
    boundary_residual = max(norm(left_residual, Inf), norm(right_residual, Inf))

    # These locations differ from MIRK's own defect-control sample points.
    defect_by_state, absolute_defect_by_state = zeros(field_count), zeros(field_count)
    for index in 1:length(solution.t)-1, fraction in (0.211324865405187, 0.788675134594813)
        position = solution.t[index] + fraction * (solution.t[index+1] - solution.t[index])
        state, derivative = solution(position), solution(position, Val{1})
        rhs = similar(state)
        callbacks.rhs!(rhs, state, nothing, position)
        difference = abs.(view(derivative, 1:field_count) - view(rhs, 1:field_count))
        defect_by_state .= max.(defect_by_state, difference ./ (1 .+ abs.(view(rhs, 1:field_count))))
        absolute_defect_by_state .= max.(absolute_defect_by_state, state_scales .* difference)
    end

    (; inner, outer) = problem.integration
    interval_count = max(1, ceil(Int, (outer - inner) / sample_step))
    positions = collect(range(inner, outer; length=interval_count+1))
    states = hcat([unpack_complex(state_scales .* view(solution(position), 1:field_count))
        for position in positions]...)
    profile = two_phase_collocation_profile(reference_offsets, amplitudes,
        problem, packed_reference, parameter_scales)
    admissible = two_phase_profile_admissible(pack_two_phase_profile(profile, problem), problem)
    finite = all(isfinite, states) && all(state -> all(isfinite, state), solution.u) &&
        all(isfinite, defect_by_state) && all(isfinite, absolute_defect_by_state) &&
        isfinite(boundary_residual) && isfinite(parameter_drift)
    passed = result.successful && finite && admissible &&
        boundary_residual <= result.settings.nonlinear_tolerance &&
        maximum(defect_by_state) <= result.settings.tolerance && parameter_drift <= parameter_tolerance
    return (; passed, profile, positions, states, boundary_residual, parameter_drift,
        defect_by_state, absolute_defect_by_state, admissible)
end
