"""
    TwoPhaseProblem{Components, Couplings, PolarValues, Grid,
        Phase1Nodes, Phase2Nodes, Integration}

Fixed physical and numerical settings shared by two-cloud profile
correction, response calculations, and evolution.

Fields:
- `spin`, `alpha`: black-hole spin and gravitational coupling.
- `components`, `cloud_m`: retained Fourier/angular components and the two
  independent clouds' azimuthal numbers.
- `couplings`, `polar_values`, `grid`: angular coupling data, polar basis
  values, and angular quadrature.
- `phase1_nodes`, `phase2_nodes`: sampling nodes for the two independent phases.
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
stores the supplied fields without performing those calculations or checks.

Trial frequencies, boundary coefficients, and current cloud amplitudes are
supplied separately; this object does not contain a solved profile.
"""
Base.@kwdef struct TwoPhaseProblem{
    Components, Couplings, PolarValues, Grid,
    Phase1Nodes, Phase2Nodes, Integration
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

    # Nonlinear potential and numerical coordinates
    potential::Symbol
    coefficient_scale::Float64
    free_outer_coefficient_indices::Vector{Int}

    # Boundary branches and radial integration settings
    initial_bound::BitVector
    integration::Integration
end

"""
    TwoPhaseProfile(frequencies, horizon_coefficients, outer_coefficients)

Frequencies and physical boundary coefficients specifying a guess at, or a
converged profile of, the two-cloud system.

- `frequencies`: the two complex cloud frequencies, in cloud order. Each
  Fourier component's frequency is their combination specified by component.q.
- `horizon_coefficients`: one physical horizon coefficient per retained radial
  mode, ordered by component and then by that component's l_values.
- `outer_coefficients`: one physical outer coefficient per retained radial
  mode, in the same order. Includes the two prescribed amplitudes at the first
  retained angular harmonic of q=(1,0) and q=(0,1).

All three fields are vectors of Complex{T}, with a shared real numeric type T.
Boundary coefficients retain their physical normalization. Use
pack_two_phase_profile(profile, problem) to select the free outer coefficients,
scale the boundary coefficients, and pack the unknowns into a real solver vector.

This object stores boundary data rather than sampled radial functions.
Construction does not integrate the radial equations or verify matching.
"""
struct TwoPhaseProfile{T<:Real}
    frequencies::Vector{Complex{T}}
    horizon_coefficients::Vector{Complex{T}}
    outer_coefficients::Vector{Complex{T}}
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
    (; spin, alpha, components, cloud_m, couplings, polar_values, grid,
        phase1_nodes, phase2_nodes, potential) = problem
    (; inner, outer, match, step, horizon_order, outer_order) = problem.integration

    inner < match < outer || throw(ArgumentError("match must lie between endpoints"))
    isempty(components) && throw(ArgumentError("require at least one component"))
    counts = [length(component.l_values) for component in components]
    offsets = cumsum(vcat(0, counts))
    length(horizon_coefficients) == length(outer_coefficients) == offsets[end] ||
        throw(DimensionMismatch("require one boundary coefficient per retained pair"))

    function assemble_boundary(side, position, coefficients)
        blocks = [component_boundary_state(side, position, spin, alpha, component,
            cloud_frequencies, cloud_m, couplings,
            coefficients[offsets[j]+1:offsets[j+1]]; horizon_order, outer_order)
            for (j, component) in enumerate(components)]

        # Each block is [R; Q]. Gather ALL R blocks before ALL Q blocks for the RHS.
        profiles = [block[1:counts[j]] for (j, block) in enumerate(blocks)]
        slopes = [block[counts[j]+1:end] for (j, block) in enumerate(blocks)]
        return vcat(profiles..., slopes...)
    end

    horizon_state = assemble_boundary(:horizon, inner, horizon_coefficients)
    outer_state = assemble_boundary(:outer, outer, outer_coefficients)
    rhs! = (derivative, state, position) -> nonlinear_two_phase_rhs!(
        derivative, state, position, spin, alpha, components, cloud_frequencies,
        cloud_m, couplings, polar_values, grid, phase1_nodes, phase2_nodes; potential)

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
