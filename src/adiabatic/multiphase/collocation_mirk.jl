"""
    two_phase_mirk_algorithm(; ad_chunksize=28, max_mesh_intervals=1200)

Return a fourth-order MIRK algorithm with a sparse, colored forward-AD Jacobian.

Use Newton corrections with directional-AD backtracking. `ad_chunksize` limits
the derivative batch width; `max_mesh_intervals` bounds adaptive mesh growth.
The optional backend loader installs the MIRK 1.19.2 library corrections before
this algorithm is used. Throw `ArgumentError` for nonpositive integer limits.
"""
function two_phase_mirk_algorithm(; ad_chunksize::Integer=28, max_mesh_intervals::Integer=1200)
    ad_chunksize > 0 && max_mesh_intervals > 0 ||
        throw(ArgumentError("require positive AD chunk and mesh limits"))
    coloring = BoundaryValueDiffEqMIRK.BoundaryValueDiffEqCore.GreedyColoringAlgorithm()
    jacobian_ad = AutoSparse(AutoForwardDiff(; chunksize=ad_chunksize);
        coloring_algorithm=coloring)
    linesearch = BackTracking(; maxiters=20, autodiff=AutoForwardDiff(; chunksize=1))
    nonlinear_solver = NewtonRaphson(; jacobian_reuse=false, linesearch)
    return MIRK4(; jac_alg=BVPJacobianAlgorithm(jacobian_ad),
        nlsolve=nonlinear_solver, max_num_subintervals=max_mesh_intervals)
end

"""
    solve_two_phase_collocation(problem::TwoPhaseProblem, amplitudes,
        reference_profile::TwoPhaseProfile, initial_guess;
        state_scales, parameter_scales, initial_step=1.0, tolerance=1e-4,
        nonlinear_tolerance=tolerance/100, max_mesh_intervals=1200,
        max_nonlinear_iterations=20, time_limit=2400.0, ad_chunksize=28, trace=false)

Solve the joint radial BVP with adaptive fourth-order MIRK collocation.

`initial_guess(position)` returns the physical complex field/slope vector
`[R; Q]` at tortoise radius `position`, throughout the problem's inner/outer
interval. Start the parameter offsets at zero relative to `reference_profile`.
Supply fixed positive `state_scales` and `parameter_scales`; the latter use
the ordering of [`pack_two_phase_profile`](@ref). Copies of these arrays and
`amplitudes` are retained for this solve.

`initial_step` sets initial mesh spacing. `tolerance` controls the scaled
differential defect; `nonlinear_tolerance` controls the collocation equations.
Limit mesh growth and Newton corrections with `max_mesh_intervals` and
`max_nonlinear_iterations`. `time_limit` bounds elapsed solve time through RHS
checks, including compilation within the solve; `trace` enables Newton output.

Return a named tuple containing the raw `solution`, decoded boundary `profile`,
`successful` library-return flag, fixed coordinate arrays, and `settings`.
Keep the associated `problem` to interpret mode ordering. Library success
requires separate physical boundary, off-mesh defect, and refinement checks
before accepting the profile. Solver exceptions propagate to the caller.
"""
function solve_two_phase_collocation(problem::TwoPhaseProblem, amplitudes,
        reference_profile::TwoPhaseProfile, initial_guess;
        state_scales, parameter_scales, initial_step=1.0, tolerance=1e-4,
        nonlinear_tolerance=tolerance/100, max_mesh_intervals::Integer=1200,
        max_nonlinear_iterations::Integer=20, time_limit=2400.0,
        ad_chunksize::Integer=28, trace::Bool=false)
    for value in (initial_step, tolerance, nonlinear_tolerance, time_limit)
        isfinite(value) && value > 0 ||
            throw(ArgumentError("require positive finite spacings, tolerances, and time limit"))
    end
    max_nonlinear_iterations > 0 || throw(ArgumentError("require a positive Newton limit"))
    (; inner, outer) = problem.integration
    isfinite(inner) && isfinite(outer) && inner < outer ||
        throw(ArgumentError("require finite ordered radial endpoints"))
    amplitudes, state_scales, parameter_scales = copy(amplitudes), copy(state_scales), copy(parameter_scales)
    packed_reference = pack_two_phase_profile(reference_profile, problem)
    coordinates = (; state_scales, packed_reference, parameter_scales)
    callbacks = two_phase_collocation_callbacks(problem, amplitudes; coordinates...)
    field_count = length(state_scales)

    function guess(unused, position)
        physical_state = initial_guess(position)
        length(physical_state) == field_count ÷ 2 ||
            throw(DimensionMismatch("initial_guess must return physical [R; Q]"))
        all(isfinite, physical_state) || throw(ArgumentError("initial_guess must be finite"))
        return vcat(pack_complex(physical_state) ./ state_scales,
            zeros(eltype(packed_reference), length(packed_reference)))
    end
    started = time()
    function rhs!(derivative, state, unused, position)
        time() - started <= time_limit || error("Radial collocation time limit reached")
        callbacks.rhs!(derivative, state, unused, position)
        return nothing
    end
    bvp = TwoPointBVProblem(rhs!, (callbacks.horizon!, callbacks.outer!), guess,
        (inner, outer); bcresid_prototype=(zeros(field_count), zeros(field_count)))
    algorithm = two_phase_mirk_algorithm(; ad_chunksize, max_mesh_intervals)
    solution = solve(bvp, algorithm; dt=initial_step, abstol=tolerance,
        nlsolve_kwargs=(; abstol=nonlinear_tolerance, reltol=0.0,
            maxiters=max_nonlinear_iterations, maxtime=time_limit, show_trace=Val(trace)))
    parameter_offsets = view(first(solution.u), (field_count + 1):length(first(solution.u)))
    profile = two_phase_collocation_profile(parameter_offsets, amplitudes,
        problem, packed_reference, parameter_scales)
    settings = (; initial_step, tolerance, nonlinear_tolerance, max_mesh_intervals,
        max_nonlinear_iterations, time_limit, ad_chunksize, trace)
    return (; solution, profile, successful=successful_retcode(solution), amplitudes,
        state_scales, parameter_scales, packed_reference, settings)
end
