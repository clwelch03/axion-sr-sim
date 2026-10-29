"""
    apply_two_phase_mirk_corrections!()

Install the tested Jacobian-band and mesh-transfer corrections for MIRK 1.19.2.

Correct the two-point Jacobian bandwidth for fully coupled states, and preserve
old nodal values while interpolating a new mesh. These method overrides apply
throughout the current Julia process and support serial BVP solves. Repeated
calls reuse the installed corrections. Return `nothing`; reject other MIRK
versions until their Jacobian and mesh-transfer behavior have been verified.
"""
function apply_two_phase_mirk_corrections!()
    Base.pkgversion(BoundaryValueDiffEqMIRK) == v"1.19.2" ||
        error("Verify the radial BVP corrections before changing MIRK versions")
    isdefined(BoundaryValueDiffEqMIRK, :_axion_radial_mirk_corrected) && return nothing

    @eval BoundaryValueDiffEqMIRK begin
        function __generate_sparse_jacobian_prototype(
                ::MIRKCache, ::TwoPointBVProblem, left, right, state_count, node_count)
            fast_scalar_indexing(left) || error("Scalar indexing required")
            left_count, right_count = length(left), length(right)
            rows = left_count + right_count + state_count * (node_count - 1)
            columns = state_count * node_count
            lower = min(rows - 1, max(left_count + state_count - 1,
                left_count + right_count - 1))
            upper = min(columns - 1, max(2state_count - left_count - 1, state_count - 1))
            prototype = BandedMatrix(Ones{eltype(left)}(rows, columns), (lower, upper))
            return rows < columns ? sparse(prototype) : prototype
        end

        const _axion_radial_old_mesh = Ref{Any}(nothing)
        const _axion_radial_old_values = Ref{Any}(nothing)
        function interp_eval!(output::AbstractArray, cache::MIRKCache, position, mesh, mesh_steps)
            if _axion_radial_old_mesh[] !== mesh
                _axion_radial_old_mesh[] = mesh
                _axion_radial_old_values[] = copy.(cache.y₀.u)
            end
            index = interval(mesh, position)
            step = mesh_steps[index]
            weights, _ = interp_weights((position - mesh[index]) / step, cache.alg)
            overwritten_base = copy(cache.y₀.u[index])
            sum_stages!(output, cache, weights, index, step)
            output .+= _axion_radial_old_values[][index] .- overwritten_base
            return output
        end
        const _axion_radial_mirk_corrected = true
    end
    return nothing
end
