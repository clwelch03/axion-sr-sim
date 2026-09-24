using Test, Random
include("model.jl")
include("correction.jl")
include("spatial.jl")
root = normpath(joinpath(@__DIR__, "../../.."))
saved = deserialize(joinpath(root, "runs/adiabatic/radialsolvertest/seed_cosine.jls"))
context = benchmark_context(saved)
solution = last(saved.continuation).solution
Random.seed!(42)

@testset "Generic physics and derivatives" begin
    for side in (:horizon, :outer)
        amplitudes = side === :horizon ? solution.horizon_amplitudes : solution.outer_amplitudes
        position = side === :horizon ? context.options.inner : context.options.outer
        @test generic_boundary(side, solution.frequency, amplitudes, context) ≈
            boundary_state(side, position, context.spin, context.alpha, solution.frequency,
                context.modes, context.couplings, amplitudes)
    end
    for position in (-50.0, 5.0, 30.0), magnitude in (1e-4, 0.2)
        state = magnitude*randn(ComplexF64, 12)
        reference = similar(state)
        nonlinear_tortoise_rhs!(reference, state, position, context.spin, context.alpha,
            solution.frequency, context.modes, context.couplings, context.polar,
            context.azimuthal, context.grid)
        geometry = tortoise_geometry(position, context.spin)
        @test radial_rhs(pack_complex(state), solution.frequency, geometry, context) ≈ pack_complex(reference)
        packed = vcat(pack_complex(state), real(solution.frequency), imag(solution.frequency))
        function_rhs = value -> radial_rhs(value[1:24], complex(value[25],value[26]), geometry, context)
        direction = randn(26)
        jacobian = ForwardDiff.jacobian(function_rhs, packed)
        centered = (function_rhs(packed+1e-6direction)-function_rhs(packed-1e-6direction))/2e-6
        @test norm(jacobian*direction-centered)/norm(centered) < 1e-7
    end
end

@testset "Correction convergence" begin
    residual = state -> [state[1]^2-2, state[2]-3]
    for strategy in (:residual, :deuflhard)
        result = correct_system(residual, state -> ForwardDiff.jacobian(residual,state),
            [10.0, 0.0]; strategy)
        @test result.status == "converged"
        @test result.state ≈ [sqrt(2), 3]
    end
end

@testset "Sparse spatial Jacobians" begin
    short_context = merge(context, (; options=(; inner=0.0, match=0.2, outer=0.4, step=0.1)))
    amplitude = last(saved.continuation).amplitude
    problem = shooting_problem(short_context, amplitude, solution)
    for (method, resolution, align_mesh) in ((:multiple, 2, false), (:multiple, 2, true), (:collocation, 0.1, false))
        spatial = spatial_problem(problem; method, resolution, align_mesh)
        matrix = spatial.jacobian(spatial.initial)
        direction = randn(length(spatial.initial))
        direction .*= max.(abs.(spatial.initial), 1e-8)
        centered = (spatial.residual(spatial.initial+1e-5direction) -
            spatial.residual(spatial.initial-1e-5direction))/2e-5
        @test size(matrix) == (length(spatial.initial), length(spatial.initial))
        @test norm(matrix*direction-centered)/norm(centered) < 1e-5
    end
end
