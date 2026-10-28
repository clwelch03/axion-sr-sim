"""Experiment-local correction for the fully coupled two-point Jacobian band.

The installed MIRK 1.19.2 assumes bandwidth M+1 on either side. With La left
boundary equations, a dense two-node collocation block instead needs lower
bandwidth La+M-1 and upper bandwidth 2M-La-1. This changes storage/coloring,
not the collocation equations. test_mirk_jacobian.jl compares both Jacobians
against uncompressed AD before this patch is used for a radial solve.
"""
function apply_mirk_coupling_patch!()
    Base.pkgversion(BoundaryValueDiffEqMIRK) == v"1.19.2" ||
        error("Recheck the Jacobian correction before changing MIRK versions")
    @eval BoundaryValueDiffEqMIRK begin
        function __generate_sparse_jacobian_prototype(
                ::MIRKCache, ::TwoPointBVProblem, ya, yb, M, N)
            fast_scalar_indexing(ya) || error("Scalar indexing required")
            La,Lb = length(ya),length(yb)
            rows,columns = La+Lb+M*(N-1), M*N
            lower = min(rows-1,max(La+M-1,La+Lb-1))
            upper = min(columns-1,max(2M-La-1,M-1))
            J = BandedMatrix(Ones{eltype(ya)}(rows,columns),(lower,upper))
            rows < columns ? sparse(J) : J
        end
    end
    nothing
end

# ADTypes.AutoSparse defaults to NoColoringAlgorithm when constructed by the
# caller. Specify coloring explicitly when overriding the AD chunk size.
mirk_sparse_ad(chunksize) = AutoSparse(AutoForwardDiff(;chunksize),
    coloring_algorithm=BoundaryValueDiffEqMIRK.BoundaryValueDiffEqCore.GreedyColoringAlgorithm())
