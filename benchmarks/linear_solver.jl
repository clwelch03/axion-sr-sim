# Compare the block factorisation with the sparse LU (UMFPACK) reference on the benchmark states.
#
#   julia --project=. benchmarks/linear_solver.jl [state names...]
#
# For each state, at the converged unknowns and at the fixed displaced unknowns (`displaced_unknowns`): times of the
# pieces of the factorisation (minimum of 5 calls), the backward error and the difference between the two solvers on a
# fixed vector, a fixed 4-column matrix and a Newton right-hand side, and a lower bound on the condition number.
# Needs the benchmark harness `common.jl` (its directory is BENCH_COMMON_DIR, default this directory).
include(joinpath(get(ENV, "BENCH_COMMON_DIR", @__DIR__), "common.jl"))
using .BenchmarkHarness
using .BenchmarkHarness: load_state, spec_named, jacobian_at, displaced_unknowns, newton_scaling, residual_at, fixed_vector
using AxionCloud, SparseArrays, LinearAlgebra, Printf
const AC = AxionCloud

best(piece; repeats=5) = (piece(); minimum(@elapsed(piece()) for _ in 1:repeats))

"""Backward error `‖M x − b‖∞ / (‖M‖∞ ‖x‖∞ + ‖b‖∞)` of an approximate solution `x` of `M x = b`."""
backward_error(matrix, x, b) = norm(matrix * x - b, Inf) / (opnorm(matrix, Inf) * norm(x, Inf) + norm(b, Inf))

"""Lower bound on `‖M‖₁ ‖M⁻¹‖₂` from a few power iterations with `solve`, which applies `M⁻¹` (no transposed solves)."""
function condition_lower_bound(solve, matrix; iterations=3)
    vector = fixed_vector(size(matrix, 1); frequency=0.91)
    growth = 0.0
    for _ in 1:iterations
        solved = solve(vector)
        growth = max(growth, norm(solved) / norm(vector))
        vector = solved ./ norm(solved)
    end
    return growth * opnorm(matrix, 1)
end

function report_state(name)
    state = load_state(spec_named(name))
    scaling = newton_scaling(state)
    layout = scaling.layout
    println("== $name: $(length(state.solution.unknowns)) unknowns, BLAS threads $(BLAS.get_num_threads())")
    for (label, unknowns) in (("converged", state.solution.unknowns), ("displaced", displaced_unknowns(state.solution)))
        jacobian = jacobian_at(state, unknowns)
        scaled = Diagonal(scaling.row) * jacobian * Diagonal(scaling.column)
        reference = AC.ScaledFactorization(jacobian, scaling; method=:sparse_lu)
        block = AC.ScaledFactorization(jacobian, scaling; method=:block)
        umfpack = lu(scaled)
        println("-- $label unknowns")
        @printf("   factorisation [s]: scaled product %.4f | UMFPACK lu %.4f, lu! reusing symbolic %.4f | block total %.4f, of which extraction %.4f\n",
                best(() -> Diagonal(scaling.row) * jacobian * Diagonal(scaling.column)), best(() -> lu(scaled)), best(() -> lu!(umfpack, scaled)),
                best(() -> AC.ScaledFactorization(jacobian, scaling; method=:block)),
                best(() -> AC.stamp_panels(jacobian, layout, scaling.row, scaling.column)))
        vector, matrix = fixed_vector(size(jacobian, 1)), hcat((fixed_vector(size(jacobian, 1); frequency=0.2k) for k in 1:4)...)
        @printf("   solve [s], vector: UMFPACK %.5f block %.5f; 4 columns: UMFPACK %.5f block %.5f\n",
                best(() -> reference \ vector), best(() -> block \ vector), best(() -> reference \ matrix), best(() -> block \ matrix))
        for (rhs_label, rhs) in (("vector", vector), ("4 columns", matrix), ("Newton rhs", -residual_at(state, unknowns)))
            solution_reference, solution_block = reference \ rhs, block \ rhs
            scaled_rhs = scaling.row .* rhs
            worst_error(solution) = maximum(backward_error(scaled, solution[:, j] ./ scaling.column, scaled_rhs[:, j]) for j in 1:size(rhs, 2))
            @printf("   %-10s backward error: UMFPACK %.2e block %.2e; |block - UMFPACK| / |UMFPACK| in scaled unknowns %.2e\n", rhs_label,
                    worst_error(solution_reference), worst_error(solution_block),
                    norm((solution_block - solution_reference) ./ scaling.column) / norm(solution_reference ./ scaling.column))
        end
        @printf("   condition number lower bound: UMFPACK %.2e block %.2e\n", condition_lower_bound(vector -> umfpack \ vector, scaled),
                condition_lower_bound(vector -> block.factorization \ vector, scaled))
    end
end

foreach(report_state, isempty(ARGS) ? ["linear", "plateau", "nonlinear", "large_mesh"] : ARGS)
