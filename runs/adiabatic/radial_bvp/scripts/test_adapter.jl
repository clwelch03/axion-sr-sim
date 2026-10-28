using Test
include("inputs.jl")

@testset "Frozen radial equations exposed as a real BVP" begin
    case = initial_case("weak")
    callbacks = collocation_callbacks(case)
    p = zeros(56)
    seed = reference_record("ad_seed_20261027_v1-seed")
    @test case.profile.frequencies == seed.profile.frequencies
    @test case.profile.outer_coefficients[[1,7]] == case.amplitudes
    for (side,branch,bc!) in ((:horizon,seed.radial_result.from_horizon,callbacks.bca!),
            (:outer,seed.radial_result.from_outer,callbacks.bcb!))
        endpoint = endpoint_state(case.profile,case.problem,side)
        @test endpoint ≈ pack_complex(branch.states[:,1]) rtol=2e-14
        residual = zeros(56)
        bc!(residual,endpoint./case.state_scales,p)
        @test norm(residual,Inf) <= 2e-14
    end
    for position in (-37.234, 41.234, 217.234)
        u = callbacks.guess(p,position)
        rhs = (state,parameters) -> begin
            buffer = similar(state, promote_type(eltype(state),eltype(parameters)))
            callbacks.rhs!(buffer,state,parameters,position)
            buffer
        end
        # The existing Dual specialization must also work when both profile
        # values and frequency/boundary parameters are BVP unknowns.
        z = vcat(u,p)
        f = z -> rhs(z[1:56],z[57:end])
        direction = sin.(1.0:112.0)
        directional_derivative = ForwardDiff.jacobian(f,z)*direction
        # Sweep steps: cancellation in the very small driven channels limits
        # finite differences before it limits the local analytic derivative.
        errors = map((1e-3,3e-4,1e-4,3e-5)) do h
            finite_difference = (f(z+h*direction)-f(z-h*direction))/(2h)
            norm(directional_derivative-finite_difference)/norm(directional_derivative)
        end
        println("directional derivative at $position: ", errors)
        @test minimum(errors) < 2e-7
        r,q = sample_block(seed.radial_result,position,0,14)
        rplus = first(sample_block(seed.radial_result,position+1e-5,0,14))
        rminus = first(sample_block(seed.radial_result,position-1e-5,0,14))
        @test norm((rplus-rminus)/2e-5-q)/norm(q) < 1e-7
    end
    augmented = vcat(callbacks.guess(p,-100.0), p)
    boundary_jacobian = ForwardDiff.jacobian(augmented) do state
        result = similar(state,56)
        callbacks.bca!(result,state,p)
        result
    end
    @test norm(boundary_jacobian[:,57:end]) > 0
    for name in ("target_zero","target_driven","scaled_weak")
        target = initial_case(name)
        @test target.amplitudes == [0.036390480353373555,1.0896058341192524e-4]
        @test all(isfinite,target.state_scales) && all(>(0),target.state_scales)
        @test all(isfinite,target.parameter_scales) && all(>(0),target.parameter_scales)
        @test all(isfinite,target.guess(40.0))
    end
end
