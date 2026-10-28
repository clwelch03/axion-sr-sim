using Serialization: deserialize, serialize

"""
    read_branch(checkpoint_path; checkpoint_adapter=nothing)

Read a ContinuationBranch stored directly in a file, or in the `branch`
field of a checkpoint named tuple that also holds restart information.

For an older format, supply
`checkpoint_adapter(saved_checkpoint; checkpoint_path)`.
That function must convert the entire checkpoint to a ContinuationBranch.

Return only the branch. This does not build a runtime problem, reconstruct
radial samples, or check convergence.
"""
function read_branch(checkpoint_path; checkpoint_adapter=nothing)
    saved_checkpoint = deserialize(checkpoint_path)

    # A file may contain just the shared branch object.
    saved_checkpoint isa ContinuationBranch && return saved_checkpoint

    branch = if saved_checkpoint isa NamedTuple && haskey(saved_checkpoint, :branch)
        # A run checkpoint can keep restart information alongside the branch.
        saved_checkpoint.branch
    elseif checkpoint_adapter !== nothing
        # The caller selects the conversion appropriate to the old file format.
        checkpoint_adapter(saved_checkpoint; checkpoint_path=abspath(checkpoint_path))
    else
        throw(ArgumentError(
            "Checkpoint does not contain a ContinuationBranch; supply a checkpoint_adapter.",
        ))
    end

    branch isa ContinuationBranch || throw(ArgumentError(
        "The checkpoint or its adapter must provide a ContinuationBranch.",
    ))
    return branch
end

"""
    write_continuation_checkpoint(path, checkpoint; overwrite=true)

Serialize to `path * ".tmp"`, then rename the completed file to `path`.
With `overwrite=false`, refuse existing output or temporary files. Preserve
the supplied checkpoint schema so existing seeds and recovery tools still work.
"""
function write_continuation_checkpoint(path, checkpoint; overwrite=true)
    temporary = path * ".tmp"
    if !overwrite && (ispath(path) || ispath(temporary))
        error("Output files already exist; choose a fresh output prefix: $path")
    end
    serialize(temporary, checkpoint)
    mv(temporary, path; force=overwrite)
    return path
end

"""
    continuation_newton_checkpoint(path, metadata; log_io=stdout)

Return a Newton `accepted_state_callback` that saves the latest accepted iterate
alongside its fixed problem and normalization metadata. This checkpoint is for
recovery; an accepted Newton iterate need not be a converged branch point.
"""
function continuation_newton_checkpoint(path, metadata; log_io=stdout)
    return function(iteration, state, mismatch)
        write_continuation_checkpoint(path, (; metadata..., iteration, state, mismatch))
        println(log_io, "Newton ", iteration, " weighted residual=", norm(mismatch))
        flush(log_io)
    end
end
