"""
    TwoPhaseProblem{Components, Couplings, PolarValues, Grid,
        Phase1Nodes, Phase2Nodes, SourceBasis, Integration}

Fixed physical and numerical settings shared by two-cloud profile
correction, response calculations, and evolution.

Fields:
- `spin`, `alpha`: black-hole spin and gravitational coupling.
- `components`, `cloud_m`: retained Fourier/angular components and the two
  independent clouds' azimuthal numbers.
- `couplings`, `polar_values`, `grid`: angular coupling data, polar basis
  values, and angular quadrature.
- `phase1_nodes`, `phase2_nodes`: sampling nodes for the two independent phases.
- `source_basis`: fixed reconstruction/projection data, built once from the
  retained mode order and angular/phase grids when this problem is constructed.
- `potential`: nonlinear potential choice.
- `coefficient_scale`: positive scale used to normalize boundary coefficients
  and matching residuals; held fixed as cloud amplitudes vary.
- `free_outer_coefficient_indices`: positions of the outer coefficients solved
  for, in the flattened component/angular-mode ordering.
- `initial_bound`: one Boolean per component, recording whether its reference
  frequency lies below alpha. Trial frequencies must preserve these branches.
- `integration`: radial step, endpoints, matching location, and boundary-series
  orders.

The driver calculates the coefficient indices and reference boundary branches,
constructs this object with keywords, and validates the settings. The constructor
stores the supplied fields and builds the source basis by default, without
performing those driver calculations or checks. A supplied basis must match
the retained mode order, cloud m values, polar table, and phase grids. Rebuild
the problem/basis if those settings change; treat its arrays as read-only.

Trial frequencies, boundary coefficients, and current cloud amplitudes are
supplied separately; this object does not contain a solved profile.
"""
Base.@kwdef struct TwoPhaseProblem{
    Components, Couplings, PolarValues, Grid,
    Phase1Nodes, Phase2Nodes, SourceBasis<:TwoPhaseSourceBasis, Integration
}
    # Background and retained modes
    spin::Float64
    alpha::Float64
    components::Components
    cloud_m::Vector{Int}

    # Angular projection and phase sampling
    couplings::Couplings
    polar_values::PolarValues
    grid::Grid
    phase1_nodes::Phase1Nodes
    phase2_nodes::Phase2Nodes
    source_basis::SourceBasis = build_two_phase_source_basis(
        [(l, component.q) for component in components for l in component.l_values],
        cloud_m, polar_values, grid, phase1_nodes, phase2_nodes)

    # Nonlinear potential and numerical coordinates
    potential::Symbol
    coefficient_scale::Float64
    free_outer_coefficient_indices::Vector{Int}

    # Boundary branches and radial integration settings
    initial_bound::BitVector
    integration::Integration
end

"""
    two_phase_coefficient_indices(components)

Return prescribed outer coefficient positions in cloud order,
followed by free positions in flattened component/angular order.
Require both fundamental components.
"""
function two_phase_coefficient_indices(components)
    labels = [(l, component.q)
        for component in components for l in component.l_values]

    fixed_outer_coefficient_indices = [findfirst(pair -> pair[2] == q, labels)
        for q in ((1, 0), (0, 1))]
    any(isnothing, fixed_outer_coefficient_indices) &&
        throw(ArgumentError("require both fundamental components"))

    fixed_outer_coefficient_indices = Int.(fixed_outer_coefficient_indices)
    free_outer_coefficient_indices = setdiff(1:length(labels), fixed_outer_coefficient_indices)
    return (; fixed_outer_coefficient_indices, free_outer_coefficient_indices)
end

"""
    pack_two_phase_profile(profile::TwoPhaseProfile, problem::TwoPhaseProblem)

Scale the boundary coefficients, prepend the complex frequencies, and return a real representation
of the complex values. (The prescribed outer amplitudes are omitted because we don't solve for them.)
"""
function pack_two_phase_profile(profile::TwoPhaseProfile, problem::TwoPhaseProblem)
    scaled_coefficients = vcat(profile.horizon_coefficients,
        profile.outer_coefficients[problem.free_outer_coefficient_indices]) ./ problem.coefficient_scale
    return pack_complex(vcat(profile.frequencies, scaled_coefficients))
end

"""
    predict_two_phase_amplitude(seed, amplitude; cloud_index=2)

Change one prescribed amplitude using weak-cloud scaling: boundary coefficients
scale by the amplitude ratio to `abs(q[cloud_index])`, and frequencies stay fixed.
Residual weights use ONE power of the ratio for every dependent mode, including
higher harmonics, as in the seed normalization. Leave the seed arrays unchanged.
"""
function predict_two_phase_amplitude(seed, amplitude; cloud_index=2)
    cloud_index in (1, 2) || throw(ArgumentError("cloud_index must be 1 or 2"))
    isfinite(amplitude) && amplitude > 0 ||
        throw(ArgumentError("require a positive finite amplitude"))
    ratio = amplitude / seed.cloud_amplitudes[cloud_index]
    mode_labels = seed.problem.source_basis.mode_labels
    factors = [ratio^abs(q[cloud_index]) for (_, q) in mode_labels]
    profile = TwoPhaseProfile(copy(seed.profile.frequencies),
        factors .* seed.profile.horizon_coefficients,
        factors .* seed.profile.outer_coefficients)
    cloud_amplitudes = copy(seed.cloud_amplitudes)
    cloud_amplitudes[cloud_index] = amplitude
    residual_weights = seed.residual_weights .* repeat(
        [q[cloud_index] == 0 ? 1.0 : ratio for (_, q) in mode_labels], 4)
    return (; profile, cloud_amplitudes, residual_weights)
end

# Physics-specific correction; stepping and persistence live in core/continuation*.
function correct_two_phase_amplitude(problem::TwoPhaseProblem, prediction;
        tolerance=1e-5, max_iterations=20, accepted_state_callback=nothing)
    (; profile, cloud_amplitudes, residual_weights) = prediction
    initial_state = pack_two_phase_profile(profile, problem)
    residual = state -> two_phase_real_matching_residual(
        state, cloud_amplitudes, problem) ./ residual_weights
    admissible = state -> two_phase_profile_admissible(state, problem)
    solution = newton_solve(residual, initial_state; tolerance, max_iterations,
        admissible, accepted_state_callback, jacobian_backend=:forwarddiff,
        jacobian_chunk_size=length(initial_state))
    profile = unpack_two_phase_profile(solution.state, cloud_amplitudes, problem)
    return (; problem, cloud_amplitudes, residual_weights, solution, profile)
end

# Run outside the retry block: a failed fresh integration must stop the run.
function validate_two_phase_amplitude(point; tolerance)
    radial_result = two_phase_matching_residual(
        point.profile, point.problem; save_radial_result=true)
    mismatch = pack_complex(radial_result.residual) ./
        point.problem.coefficient_scale ./ point.residual_weights
    weighted_residual_norm = norm(mismatch)
    weighted_residual_norm <= tolerance || error(
        "Fresh integration exceeds the matching tolerance: $weighted_residual_norm")
    all(isfinite, radial_result.from_horizon.states) &&
        all(isfinite, radial_result.from_outer.states) ||
        error("Fresh integration contains nonfinite radial samples")
    return (; point..., radial_result, weighted_residual_norm)
end

"""
    continue_two_phase_amplitude(seed, output_prefix; target_amplitude, kwargs...)

Continue one cloud's amplitude from a saved matched seed, keeping the other
amplitude and the problem fixed. The seed contains `problem`, `cloud_amplitudes`,
`residual_weights`, `profile`, and `weighted_residual_norm`.

Each attempt saves its latest Newton iterate to `output_prefix-attemptN-accepted.jls`.
Only a converged, freshly reintegrated point is saved to `output_prefix-attemptN.jls`,
using the original seed schema. Refuse existing attempt files. To restart, supply
a matched point and a fresh prefix; return the last matched point.

The core driver controls adaptive steps. `cloud_index=2`, `tolerance=1e-5`, and
`max_iterations=20` preserve the two-cloud campaign defaults and full-width AD.
Other keywords configure [`continue_amplitude`](@ref).
"""
function continue_two_phase_amplitude(seed, output_prefix; target_amplitude,
        cloud_index=2, tolerance=1e-5, max_iterations=20, log_io=stdout, kwargs...)
    cloud_index in (1, 2) || throw(ArgumentError("cloud_index must be 1 or 2"))
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("require a positive finite tolerance"))
    check_cloud_amplitudes_and_coefficient_scale_validity(
        seed.cloud_amplitudes, seed.problem.coefficient_scale)
    isfinite(seed.weighted_residual_norm) && seed.weighted_residual_norm <= tolerance ||
        throw(ArgumentError("seed must satisfy the weighted matching tolerance"))
    length(seed.residual_weights) == 4length(seed.problem.source_basis.mode_labels) ||
        throw(DimensionMismatch("require four residual weights per radial mode"))
    all(value -> isfinite(value) && value > 0, seed.residual_weights) ||
        throw(ArgumentError("require positive finite residual weights"))

    attempt_prefix(attempt) = output_prefix * "-attempt$(attempt)"
    function correct(previous, target, attempt)
        prefix = attempt_prefix(attempt)
        any(ispath(prefix * suffix) for suffix in
            ("-accepted.jls", "-accepted.jls.tmp", ".jls", ".jls.tmp")) &&
            error("Output files already exist; choose a fresh output prefix: $prefix")
        prediction = predict_two_phase_amplitude(previous, target; cloud_index)
        metadata = (; problem=seed.problem,
            prediction.cloud_amplitudes, prediction.residual_weights)
        callback = continuation_newton_checkpoint(
            prefix * "-accepted.jls", metadata; log_io)
        return correct_two_phase_amplitude(seed.problem, prediction;
            tolerance, max_iterations, accepted_state_callback=callback)
    end
    function save(point, attempt)
        path = write_continuation_checkpoint(attempt_prefix(attempt) * ".jls",
            point; overwrite=false)
        println(log_io, "Saved matched point: ", path,
            " weighted residual=", point.weighted_residual_norm)
        flush(log_io)
    end
    return continue_amplitude(seed;
        amplitude=point -> point.cloud_amplitudes[cloud_index], correct, save,
        validate=point -> validate_two_phase_amplitude(point; tolerance),
        target_amplitude, log_io, kwargs...)
end

"""
    two_phase_boundary_state(profile::TwoPhaseProfile,
        problem::TwoPhaseProblem, side::Symbol)

Return the complex boundary-series state `[R; Q]`, where `Q = dR/dr*`.

Choose `side = :horizon` or `:outer`. Read the frequencies and the selected
boundary coefficients from `profile`; read the endpoint radius, series orders,
component ordering, and physical settings from `problem`.

Collect all field values before all slopes, ordered by `problem.components`
and then each component's `l_values`. Return two complex entries per retained
radial mode. Input arrays are read without mutation.

Throw `ArgumentError` for an unsupported side or an empty component list,
and `DimensionMismatch` if the selected coefficient count is inconsistent.
"""
function two_phase_boundary_state(profile::TwoPhaseProfile,
        problem::TwoPhaseProblem, side::Symbol)
    side in (:horizon, :outer) || throw(ArgumentError("unknown boundary side"))
    (; spin, alpha, components, cloud_m, couplings) = problem
    (; inner, outer, horizon_order, outer_order) = problem.integration
    position = side == :horizon ? inner : outer
    coefficients = side == :horizon ?
        profile.horizon_coefficients : profile.outer_coefficients

    isempty(components) && throw(ArgumentError("require at least one component"))
    counts = [length(component.l_values) for component in components]
    offsets = cumsum(vcat(0, counts))
    length(coefficients) == offsets[end] ||
        throw(DimensionMismatch("require one boundary coefficient per retained pair"))

    blocks = [component_boundary_state(side, position, spin, alpha, component,
        profile.frequencies, cloud_m, couplings,
        coefficients[offsets[index]+1:offsets[index+1]];
        horizon_order, outer_order)
        for (index, component) in enumerate(components)]

    # Each component returns [R; Q]; the joint RHS needs every R before every Q.
    fields = [block[1:counts[index]] for (index, block) in enumerate(blocks)]
    slopes = [block[counts[index]+1:end] for (index, block) in enumerate(blocks)]
    return vcat(fields..., slopes...)
end

"""
    two_phase_matching_residual(profile::TwoPhaseProfile,
        problem::TwoPhaseProblem; save_radial_result=false)

Integrate from both boundaries and return their mismatch at `match`.

`profile` supplies the complex cloud frequencies and physical boundary
coefficients. Supply one coefficient per retained (l, q), ordered by
`problem.components` and then each component's l_values.
`problem` supplies the background, angular data, phase nodes, and potential.
This function evaluates a trial; it does not fix normalizations or run Newton.

`problem.integration` must contain `inner`, `outer`, `match`, `step`,
`horizon_order`, and `outer_order`. The endpoints and matching position are
tortoise radii; require inner < match < outer. `step` is the maximum radial
RK4 spacing. Return [R_inner-R_outer; Q_inner-Q_outer] in physical normalization.
With `save_radial_result=true`, return (from_horizon, from_outer, residual), keeping
both radial integration histories and their separate matching-point values.
"""
function two_phase_matching_residual(profile::TwoPhaseProfile, problem::TwoPhaseProblem; save_radial_result=false)
    #TODO: "make the boundary-data cleanup as a separate, small refactor" - see Review recent refactoring work
    cloud_frequencies = profile.frequencies
    (; horizon_coefficients, outer_coefficients) = profile
    (; spin, alpha, components, cloud_m, couplings, source_basis, potential) = problem
    (; inner, outer, match, step, horizon_order, outer_order) = problem.integration

    inner < match < outer || throw(ArgumentError("match must lie between endpoints"))
    isempty(components) && throw(ArgumentError("require at least one component"))
    counts = [length(component.l_values) for component in components]
    offsets = cumsum(vcat(0, counts))
    length(horizon_coefficients) == length(outer_coefficients) == offsets[end] ||
        throw(DimensionMismatch("require one boundary coefficient per retained pair"))

    horizon_state = two_phase_boundary_state(profile, problem, :horizon)
    outer_state = two_phase_boundary_state(profile, problem, :outer)

    rhs! = (derivative, state, position) -> nonlinear_two_phase_rhs!(
        derivative, state, position, spin, alpha, components, cloud_frequencies,
        cloud_m, couplings, source_basis; potential)

    from_horizon = integrate_radial(rhs!, horizon_state, inner, match; step, save_profile=save_radial_result)
    from_outer = integrate_radial(rhs!, outer_state, outer, match; step, save_profile=save_radial_result)

    if save_radial_result
        residual = from_horizon.states[:, end] - from_outer.states[:, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon - from_outer
end

"""
    check_cloud_amplitudes_and_coefficient_scale_validity(cloud_amplitudes, coefficient_scale)

Require exactly two cloud amplitudes. Both amplitudes and coefficient_scale must be positive,
finite real numbers; throw otherwise.
"""
function check_cloud_amplitudes_and_coefficient_scale_validity(cloud_amplitudes, coefficient_scale)
    length(cloud_amplitudes) == 2 ||
        throw(DimensionMismatch("require two normalizations"))
    all(value -> value isa Real && isfinite(value) && value > 0, cloud_amplitudes) ||
        throw(ArgumentError("require positive finite real normalizations"))
    coefficient_scale isa Real && isfinite(coefficient_scale) && coefficient_scale > 0 ||
        throw(ArgumentError("require a positive finite real coefficient scale"))
end

"""
    unpack_two_phase_profile(packed_state, cloud_amplitudes, problem::TwoPhaseProblem)

Reconstruct a `TwoPhaseProfile` containing cloud frequencies and physical
horizon and outer boundary coefficients from a packed real candidate vector.

`packed_state` contains all real parts followed by all imaginary parts.
After `unpack_complex`, the complex coordinates contain, in order:
- Two cloud frequencies.
- All horizon boundary coefficients, divided by `problem.coefficient_scale`.
- The free outer boundary coefficients, divided by `problem.coefficient_scale`.

Coefficient ordering follows `problem.components`, then each component's
`l_values`. Multiply the stored boundary coefficients by
`problem.coefficient_scale` to recover their physical values; frequencies are
not scaled. Insert the two prescribed `cloud_amplitudes` into the outer
coefficient array at the first retained angular harmonic of q=(1,0) and q=(0,1),
respectively.

For example, with one radial mode per component, ordered as |211>, |322>,
and the driven components with frequencies 2ω1-ω2 and 2ω2-ω1, there are eight
boundary coefficients. The two fundamental outer coefficients are prescribed
as A1 and A2. Return a `TwoPhaseProfile` with fields

```text
frequencies          = [ω1, ω2]
horizon_coefficients = [h1, h2, h3, h4]
outer_coefficients   = [A1, A2, b3, b4]
```
"""
function unpack_two_phase_profile(packed_state::AbstractVector{<:Real}, cloud_amplitudes, problem::TwoPhaseProblem)
    check_cloud_amplitudes_and_coefficient_scale_validity(cloud_amplitudes, problem.coefficient_scale)
    num_radial_modes = sum(length(component.l_values) for component in problem.components)
    length(packed_state) ==  4 * num_radial_modes ||
        throw(DimensionMismatch("require two frequencies and all free coefficients")) # 4N real entries

    # Unpack the packed state into a complex vector and identify the coefficient indices
    frequency_and_boundary_guess = unpack_complex(packed_state)
    value_type = promote_type(eltype(frequency_and_boundary_guess),
                                eltype(cloud_amplitudes),
                                typeof(problem.coefficient_scale))
    (; fixed_outer_coefficient_indices, free_outer_coefficient_indices) =
        two_phase_coefficient_indices(problem.components)

    # Extract frequencies and coefficients from frequency_and_boundary_guess
    cloud_frequencies = value_type.(frequency_and_boundary_guess[1:2])
    horizon_coefficients = problem.coefficient_scale * value_type.(frequency_and_boundary_guess[3:num_radial_modes+2])
    outer_coefficients = zeros(value_type, num_radial_modes)
    outer_coefficients[fixed_outer_coefficient_indices] .= cloud_amplitudes
    outer_coefficients[free_outer_coefficient_indices] .=
        problem.coefficient_scale * frequency_and_boundary_guess[num_radial_modes+3:end]

    return TwoPhaseProfile(cloud_frequencies, horizon_coefficients, outer_coefficients)
end


"""
    solve_two_phase_mode(frequency_guess, normalizations, spin, alpha,
        components, cloud_m, couplings, polar_values, grid, phase1_nodes,
        phase2_nodes; horizon_guess, outer_guess, kwargs...)

Solve matching at two prescribed positive outer coefficients.

`normalizations[a]` fixes the first retained l in fundamental q=(1,0) or (0,1).
Supply full boundary-coefficient guess vectors in the RHS pair order.
The two fixed entries of `outer_guess` are replaced by `normalizations`.

Solve for two complex frequencies, all horizon coefficients, and free outer
coefficients. Use validated components with distinct nonzero q labels,
one per conjugate pair.
    
Scale coefficients and matching residuals by the smaller normalization;
`tolerance` applies to that scaled residual norm. Keep both fundamentals bound
and every component on its initial positive-frequency boundary branch.
Store integration keywords in the fixed `TwoPhaseProblem` used for matching.

Return frequencies, physical boundary coefficients, residuals, and iterations.
This numerical prototype does not derive a nonlinear amplitude-evolution law.
"""
function solve_two_phase_mode(frequency_guess, normalizations, spin, alpha,
        components, cloud_m, couplings, polar_values, grid, phase1_nodes, phase2_nodes;
        horizon_guess, outer_guess, potential=:cosine, tolerance=1e-9,
        max_iterations=15, difference_step=1e-6, progress_callback=nothing,
        step=0.05, inner=-100.0, outer=100.0, match=5.0,
        horizon_order=3, outer_order=7)
    length(frequency_guess) == length(normalizations) == length(cloud_m) == 2 ||
        throw(DimensionMismatch("require two frequencies, normalizations, and cloud m values"))
    all(value -> value isa Real && isfinite(value) && value > 0, normalizations) ||
        throw(ArgumentError("require positive finite real outer normalizations"))

    pairs = [(l, component.q) for component in components for l in component.l_values]
    count = length(pairs)
    length(horizon_guess) == length(outer_guess) == count ||
        throw(DimensionMismatch("require a full boundary guess at each endpoint"))
    (; free_outer_coefficient_indices) =
        two_phase_coefficient_indices(components)
    scale = minimum(normalizations)

    initial_bound = BitVector([real(sum(component.q .* frequency_guess)) < alpha
        for component in components])
    integration = (; step, inner, outer, match, horizon_order, outer_order)
    problem = TwoPhaseProblem(; spin, alpha, components, cloud_m, couplings,
        polar_values, grid, phase1_nodes, phase2_nodes, potential,
        coefficient_scale=scale, free_outer_coefficient_indices, initial_bound,
        integration)

    # Newton's vector contains frequencies, scaled horizon coefficients, then
    # scaled free outer coefficients. real_newton supplies complex candidates;
    # the profile unpacker accepts their packed real representation.
    unpack_profile = state -> unpack_two_phase_profile(
        pack_complex(state), normalizations, problem)
    initial = complex.(vcat(frequency_guess, horizon_guess ./ scale,
                        outer_guess[free_outer_coefficient_indices] ./ scale))

    function admissible(state)
        all(isfinite, state) && all(value -> 0 < real(value) < alpha, state[1:2]) || return false
        frequencies = [real(sum(component.q .* state[1:2])) for component in components]
        return all(value -> value > 0 && value != alpha, frequencies) &&
            (frequencies .< alpha) == initial_bound
    end

    residual = state -> two_phase_matching_residual(unpack_profile(state), problem) ./ scale

    # A real Jacobian is required because field reconstruction uses conjugates.
    solution = real_newton(residual, initial; admissible, tolerance, max_iterations,
        difference_step, progress_callback)

    (; frequencies, horizon_coefficients, outer_coefficients) = unpack_profile(solution.state)
    physical_residual = scale .* solution.residual
    return (; frequencies, horizon_coefficients, outer_coefficients,
        residual=physical_residual, residual_norm=norm(physical_residual),
        scaled_residual_norm=norm(solution.residual), iterations=solution.iterations)
end
