"""
    two_phase_collocation_rhs!(scaled_derivative, scaled_state, position,
        cloud_frequencies, problem::TwoPhaseProblem, state_scales)

Fill `scaled_derivative` with radial derivatives in real, scaled coordinates.

`position` is the tortoise radius in M=1 units. Supply
`scaled_state = pack_complex([R; Q]) ./ state_scales`, where `Q = dR/dr*`,
and fixed positive `state_scales` in `[Re R; Re Q; Im R; Im Q]` order.
`cloud_frequencies` supplies the complex fundamental frequencies, and `problem`
supplies the physical equations and retained-mode ordering.

Return `scaled_derivative`. Its element type must support derivatives computed
from both the state and frequencies. Supply separate state and derivative
buffers; the input state is read without mutation.
Throw `DimensionMismatch` if the state, derivative, and scale lengths differ.
"""
function two_phase_collocation_rhs!(scaled_derivative, scaled_state, position,
        cloud_frequencies, problem::TwoPhaseProblem, state_scales)
    length(scaled_derivative) == length(scaled_state) == length(state_scales) ||
        throw(DimensionMismatch("require one scale and derivative per real state entry"))

    physical_state = unpack_complex(state_scales .* scaled_state)
    derivative_type = promote_type(eltype(physical_state), eltype(cloud_frequencies))
    physical_derivative = similar(physical_state, derivative_type)

    # Calculate the derivative for the provided scaled state
    nonlinear_two_phase_rhs!(physical_derivative, physical_state, position,
        problem.spin, problem.alpha, problem.components, cloud_frequencies,
        problem.cloud_m, problem.couplings, problem.source_basis;
        potential=problem.potential)

    scaled_derivative .= pack_complex(physical_derivative) ./ state_scales
end

"""
    two_phase_collocation_profile(parameter_offsets, amplitudes,
        problem::TwoPhaseProblem, packed_reference, parameter_scales)

Return a `TwoPhaseProfile` reconstructed from scaled parameter offsets.

Supply `packed_reference = pack_two_phase_profile(reference_profile, problem)`
and fixed positive `parameter_scales` in the same packed ordering. Reconstruct
the packed parameters as `packed_reference .+ parameter_scales .* parameter_offsets`,
then unpack their frequencies and physical boundary coefficients using `problem`.
Insert the two prescribed outer coefficients from `amplitudes` in cloud order.

Zero offsets recover the reference frequencies and free coefficients. The
returned profile stores boundary data; radial field samples are supplied
separately in the collocation state. Input arrays are read without mutation.
Throw `DimensionMismatch` if the offset, reference, and scale lengths differ.
"""
function two_phase_collocation_profile(parameter_offsets, amplitudes,
        problem::TwoPhaseProblem, packed_reference, parameter_scales)
    length(parameter_offsets) == length(packed_reference) == length(parameter_scales) ||
        throw(DimensionMismatch("require one reference and scale per parameter offset"))

    packed_parameters = packed_reference .+ parameter_scales .* parameter_offsets
    return unpack_two_phase_profile(packed_parameters, amplitudes, problem)
end

"""
    two_phase_augmented_rhs!(derivative, state, position, amplitudes, problem::TwoPhaseProblem;
        state_scales, packed_reference, parameter_scales)

Fill `derivative` for the augmented radial state at tortoise radius `position`.

Supply `state = [scaled_state; parameter_offsets]`, with the field block in
`[Re R; Re Q; Im R; Im Q]` order. `state_scales` defines that block's length.
`packed_reference` and `parameter_scales` define the parameter coordinates used
by [`two_phase_collocation_profile`](@ref); `amplitudes` fixes the two prescribed
outer coefficients.

Evaluate the field derivatives using the decoded frequencies and set all
parameter derivatives to zero, enforcing radial constancy of those unknowns.
Supply separate state and derivative buffers. Return `nothing` after filling
`derivative`; read `state` without mutation. Throw `DimensionMismatch` if the
buffers do not contain the required field and parameter blocks.
"""
function two_phase_augmented_rhs!(derivative, state, position, amplitudes, problem::TwoPhaseProblem;
        state_scales, packed_reference, parameter_scales)
    field_count = length(state_scales)
    parameter_count = length(packed_reference)
    length(derivative) == length(state) == field_count + parameter_count ||
        throw(DimensionMismatch("require field entries followed by parameter offsets"))

    field_indices = 1:field_count
    parameter_indices = (field_count + 1):length(state)
    profile = two_phase_collocation_profile(view(state, parameter_indices),
        amplitudes, problem, packed_reference, parameter_scales)

    two_phase_collocation_rhs!(view(derivative, field_indices),
        view(state, field_indices), position, profile.frequencies,
        problem, state_scales)
    derivative[parameter_indices] .= 0
    return nothing
end

"""
    two_phase_collocation_boundary!(residual, state, side, amplitudes, problem::TwoPhaseProblem;
        state_scales, packed_reference, parameter_scales)

Fill `residual` with the scaled field mismatch at one radial endpoint.

Choose `side = :horizon` or `:outer`; `problem.integration` supplies the
corresponding endpoint radius and boundary-series orders. Supply
`state = [scaled_state; parameter_offsets]` and the same fixed coordinate
arrays and prescribed `amplitudes` used by [`two_phase_augmented_rhs!`](@ref).

Decode the parameter offsets, evaluate [`two_phase_boundary_state`](@ref),
and form `scaled_state - pack_complex(boundary_state) ./ state_scales`.
`residual` has one real entry per field-block entry, in
`[Re R; Re Q; Im R; Im Q]` order.

Supply separate state and residual buffers. Return `nothing` after filling
`residual`; read `state` without mutation.
Throw `DimensionMismatch` for incompatible buffer lengths and `ArgumentError`
for an unsupported boundary side.
"""
function two_phase_collocation_boundary!(residual, state, side, amplitudes, problem::TwoPhaseProblem;
        state_scales, packed_reference, parameter_scales)
    field_count = length(state_scales)
    length(residual) == field_count ||
        throw(DimensionMismatch("require one residual per field entry"))
    length(state) == field_count + length(packed_reference) ||
        throw(DimensionMismatch("require fields followed by parameter offsets"))

    parameter_offsets = view(state, (field_count + 1):length(state))
    profile = two_phase_collocation_profile(parameter_offsets, amplitudes,
        problem, packed_reference, parameter_scales)
    boundary_state = two_phase_boundary_state(profile, problem, side)
    residual .= view(state, 1:field_count) .-
        pack_complex(boundary_state) ./ state_scales
    return nothing
end

"""
    two_phase_collocation_callbacks(problem::TwoPhaseProblem, amplitudes;
        state_scales, packed_reference, parameter_scales)

Return a named tuple `(; rhs!, horizon!, outer!)` of in-place BVP callbacks.

Bind `problem`, the two prescribed outer `amplitudes`, and the fixed coordinate
arrays to these signatures:

- `rhs!(derivative, state, unused, position)` evaluates the augmented radial RHS.
- `horizon!(residual, state, unused)` evaluates the horizon endpoint residual.
- `outer!(residual, state, unused)` evaluates the outer endpoint residual.

Each callback fills its output buffer and returns `nothing`. `unused` is the
BVP library's parameter argument; the unknown parameters reside in `state`.
The captured inputs are retained by reference; keep them fixed while using
the callbacks.

For N retained radial modes, require `4N` entries in each coordinate array,
finite reference parameters, positive finite scales, and valid prescribed
amplitudes. Throw `DimensionMismatch` for incompatible coordinate or amplitude
lengths and `ArgumentError` for an empty mode set, invalid scales, nonfinite
reference values, or invalid amplitude values.
"""
function two_phase_collocation_callbacks(problem::TwoPhaseProblem, amplitudes;
        state_scales, packed_reference, parameter_scales)
    field_count = 4 * length(problem.source_basis.mode_labels)
    field_count > 0 || throw(ArgumentError("require retained radial modes"))
    length(state_scales) == length(packed_reference) == length(parameter_scales) == field_count ||
        throw(DimensionMismatch("require 4N field scales and 4N parameter coordinates"))
    for scales in (state_scales, parameter_scales)
        all(isfinite, scales) && all(>(0), scales) ||
            throw(ArgumentError("require positive finite coordinate scales"))
    end
    all(isfinite, packed_reference) ||
        throw(ArgumentError("require finite reference parameters"))
    check_cloud_amplitudes_and_coefficient_scale_validity(amplitudes, problem.coefficient_scale)

    coordinates = (; state_scales, packed_reference, parameter_scales)
    rhs! = (derivative, state, unused, position) ->
        two_phase_augmented_rhs!(derivative, state, position,
            amplitudes, problem; coordinates...)
    horizon! = (residual, state, unused) ->
        two_phase_collocation_boundary!(residual, state, :horizon,
            amplitudes, problem; coordinates...)
    outer! = (residual, state, unused) ->
        two_phase_collocation_boundary!(residual, state, :outer,
            amplitudes, problem; coordinates...)
    return (; rhs!, horizon!, outer!)
end

"""
    two_phase_collocation_scales(problem::TwoPhaseProblem, amplitudes,
        reference_profile::TwoPhaseProfile, scale_guess;
        sample_count=801, field_floor=1e-8, slope_length=100.0, frequency_scale=1e-3)

Return fixed `(; state_scales, parameter_scales)` for a collocation solve.

`scale_guess(position)` returns physical complex `[R; Q]` throughout the
problem's radial interval. Sample its per-mode peak magnitudes. `field_floor`
is a positive scalar or one positive floor per retained mode; slope scales
are at least the field scales divided by the positive length `slope_length`.
Use a populated reference for `scale_guess` when the solve starts with zero
driven fields, so both initial guesses can use the same error weights.

Normalize parameter displacements by their effect on the two scaled endpoint
states. Set the four real frequency displacement scales to `frequency_scale`.
Parameter ordering follows [`pack_two_phase_profile`](@ref). Keep the returned
scales fixed throughout Newton correction and mesh refinement.

Require at least two samples, finite ordered endpoints, positive finite floors
and scale settings, and finite guesses of the correct length. Throw
`ArgumentError` for invalid values and `DimensionMismatch` for invalid lengths.
"""
function two_phase_collocation_scales(problem::TwoPhaseProblem, amplitudes,
        reference_profile::TwoPhaseProfile, scale_guess;
        sample_count::Integer=801, field_floor=1e-8, slope_length=100.0, frequency_scale=1e-3)
    mode_count = length(problem.source_basis.mode_labels)
    floors = field_floor isa Real ? fill(field_floor, mode_count) : collect(field_floor)
    length(floors) == mode_count || throw(DimensionMismatch("require one floor per radial mode"))
    sample_count >= 2 || throw(ArgumentError("require at least two scale samples"))
    all(isfinite, floors) && all(>(0), floors) ||
        throw(ArgumentError("require positive finite field floors"))
    all(value -> isfinite(value) && value > 0, (slope_length, frequency_scale)) ||
        throw(ArgumentError("require positive finite slope and frequency scales"))
    (; inner, outer) = problem.integration
    isfinite(inner) && isfinite(outer) && inner < outer ||
        throw(ArgumentError("require finite ordered radial endpoints"))

    # Field and slope scales come from the physical reference profile.
    peaks = zeros(2mode_count)
    for position in range(inner, outer; length=sample_count)
        physical_state = scale_guess(position)
        length(physical_state) == 2mode_count ||
            throw(DimensionMismatch("scale_guess must return physical [R; Q]"))
        all(isfinite, physical_state) || throw(ArgumentError("scale_guess must be finite"))
        peaks .= max.(peaks, abs.(physical_state))
    end
    field_scales = max.(peaks[1:mode_count], floors)
    slope_scales = max.(peaks[mode_count+1:end], field_scales ./ slope_length)
    state_scales = repeat(vcat(field_scales, slope_scales), 2)

    # Parameter scales come from their effect on the scaled endpoint states.
    packed_reference = pack_two_phase_profile(reference_profile, problem)
    boundary_map = packed -> begin
        profile = unpack_two_phase_profile(packed, amplitudes, problem)
        vcat(pack_complex(two_phase_boundary_state(profile, problem, :horizon)) ./ state_scales,
            pack_complex(two_phase_boundary_state(profile, problem, :outer)) ./ state_scales)
    end
    boundary_jacobian = ForwardDiff.jacobian(boundary_map, packed_reference)
    influences = vec(sqrt.(sum(abs2, boundary_jacobian; dims=1)))
    parameter_scales = 1 ./ max.(influences, 1e-30)
    complex_count = length(packed_reference) ÷ 2
    parameter_scales[[1, 2, complex_count+1, complex_count+2]] .= frequency_scale
    return (; state_scales, parameter_scales)
end
