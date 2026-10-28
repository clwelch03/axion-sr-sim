"""Preserve old nodal values while MIRK interpolates onto its new mesh.

MIRK 1.19.2 writes the new guesses into cache.y₀ while subsequent interpolations
still read old values from that same array. Cache a copy for each remeshing
pass. This experiment runs one solve at a time within each Julia process.
"""
function apply_mirk_remesh_patch!()
    Base.pkgversion(BoundaryValueDiffEqMIRK) == v"1.19.2" ||
        error("Recheck mesh interpolation before changing MIRK versions")
    @eval BoundaryValueDiffEqMIRK begin
        const _radial_bvp_old_mesh = Ref{Any}(nothing)
        const _radial_bvp_old_values = Ref{Any}(nothing)
        function interp_eval!(z::AbstractArray,cache::MIRKCache,t,mesh,mesh_dt)
            if _radial_bvp_old_mesh[] !== mesh
                _radial_bvp_old_mesh[] = mesh
                _radial_bvp_old_values[] = copy.(cache.y₀.u)
            end
            i = interval(mesh,t)
            dt = mesh_dt[i]
            w,_ = interp_weights((t-mesh[i])/dt,cache.alg)
            overwritten_base = copy(cache.y₀.u[i])
            sum_stages!(z,cache,w,i,dt)
            z .+= _radial_bvp_old_values[][i] .- overwritten_base
            z
        end
    end
    nothing
end
