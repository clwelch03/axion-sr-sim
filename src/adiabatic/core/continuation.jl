using LinearAlgebra

"""
    ContinuationPoint(point_id, parameter_values, profile, diagnostics,
                      backgrounds, radial_result)

One accepted stationary calculation in a continuation branch.

The caller must check convergence before creating this record. Construction
stores the supplied objects without solving, validating, or copying them.
Treat the stored profile and diagnostic arrays as read-only after acceptance.
"""
struct ContinuationPoint{Parameters,Profile,Diagnostics,Backgrounds}
    # Identifies this calculation within its branch; not an angular/radial index.
    point_id::Int

    # Named imposed values, including amplitudes held fixed along the branch.
    # For our purposes, these will be the prescribed cloud amplitudes.
    parameter_values::Parameters

    # Existing physical profile: frequencies and both boundary-coefficient arrays.
    profile::Profile

    # Matching residuals, their scaling, actual tolerance, and solver diagnostics.
    diagnostics::Diagnostics

    # References identifying any saved branch points used as fixed backgrounds.
    backgrounds::Backgrounds

    # Existing (from_horizon, from_outer, residual) record, when samples were saved.
    radial_result::Union{Nothing,NamedTuple}
end

"""
    ContinuationBranch(settings, points, provenance)

An ordered collection of accepted stationary calculations with shared settings.
Each point records its own imposed parameter values and convergence diagnostics.
Construction stores the supplied objects without copying or validating them.
"""
struct ContinuationBranch{Settings,Point<:ContinuationPoint,Provenance}
    # Model, cloud/component ordering, fixed physical parameters, and grid settings.
    # Derived tables and solver workspaces belong in the runtime problem instead.
    settings::Settings

    # Accepted calculations in continuation order; point_id is not a vector index.
    points::Vector{Point}

    # Where the data came from, such as checkpoint paths and source revisions.
    provenance::Provenance
end

# Retry numerical failures only. Programming, checkpoint, and validation errors
# must reach the caller instead of silently shrinking the continuation step.
continuation_numerical_failure(exception) =
    exception isa SingularException || exception isa DomainError ||
    (exception isa ErrorException && (startswith(exception.msg, "Newton ") ||
        exception.msg == "nonfinite Newton residual"))

"""
    continue_amplitude(seed; amplitude, correct, save, target_amplitude, kwargs...)

Continue one positive amplitude, holding the model's other parameters fixed.
`amplitude(point)` reads its value; `correct(previous, target, attempt)` predicts
and solves a new point; `validate(point)` checks it and returns the point to
save; `save(point, attempt)` persists it. Callbacks must not mutate the previous
accepted point. Return the last saved point.

Numerical correction failures halve the actual attempted increment and retry
from the last accepted point. Validation and save errors propagate without retry.
Advance only after `save`
succeeds. The final increment is clipped to the target, even below the minimum
fractional step; a failed attempt below that minimum stops the continuation.

`iterations(point)` supplies the correction count: at most four grows the next
increment by 1.5, twelve or more halves it, and other counts leave it unchanged.
Attempt numbers start at one for each call; callers choose a fresh output prefix.
"""
function continue_amplitude(seed; amplitude, correct, save, target_amplitude,
        validate=identity,
        fractional_step=0.015625, minimum_fractional_step=1e-4,
        maximum_fractional_step=0.25,
        iterations=point -> point.solution.iterations, log_io=stdout)
    for (name, value) in ((:target_amplitude, target_amplitude),
            (:fractional_step, fractional_step),
            (:minimum_fractional_step, minimum_fractional_step),
            (:maximum_fractional_step, maximum_fractional_step))
        isfinite(value) && value > 0 ||
            throw(ArgumentError("require positive finite $name"))
    end
    minimum_fractional_step <= fractional_step <= maximum_fractional_step ||
        throw(ArgumentError("initial fractional step must lie within its bounds"))
    current_amplitude = amplitude(seed)
    isfinite(current_amplitude) && 0 < current_amplitude <= target_amplitude ||
        throw(ArgumentError("require a positive seed amplitude no larger than the target"))

    previous = seed
    attempt = 0
    while current_amplitude < target_amplitude
        attempt += 1
        target = min(target_amplitude, current_amplitude * (1 + fractional_step))
        used_fraction = (target - current_amplitude) / current_amplitude
        used_fraction > 0 || error("Continuation step cannot advance the amplitude")
        println(log_io, "Attempt ", attempt, ": amplitude ", current_amplitude, " -> ", target)
        flush(log_io)

        point = try
            correct(previous, target, attempt)
        catch exception
            continuation_numerical_failure(exception) || rethrow()
            fractional_step = used_fraction / 2
            println(log_io, "Correction failed; retrying with fractional step ", fractional_step)
            flush(log_io)
            fractional_step >= minimum_fractional_step || rethrow()
            continue
        end

        point = validate(point)
        amplitude(point) == target || error("Corrected point has the wrong amplitude")
        count = iterations(point)
        growth = count <= 4 ? 1.5 : count >= 12 ? 0.5 : 1.0
        next_fraction = clamp(growth * used_fraction,
            minimum_fractional_step, maximum_fractional_step)
        save(point, attempt)

        previous = point
        current_amplitude = target
        fractional_step = next_fraction
    end
    return previous
end
