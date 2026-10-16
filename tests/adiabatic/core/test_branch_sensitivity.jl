using LinearAlgebra: SingularException

@testset "Implicit derivative of an analytic nonlinear branch" begin
    residual = (state, amplitude) -> [
        state[1]^2 - amplitude,
        state[2] - amplitude * state[1],
    ]

    # Positive branch: u₁ = sqrt(A), u₂ = A^(3/2).
    for amplitude in (1.0, 4.0, 9.0)
        state = [sqrt(amplitude), amplitude * sqrt(amplitude)]
        @test residual(state, amplitude) == zeros(2)
        expected = [1 / (2sqrt(amplitude)), 3sqrt(amplitude) / 2]
        @test implicit_branch_derivative(residual, state, amplitude) ≈ expected
    end
end

@testset "Total energy derivative includes state and amplitude dependence" begin
    residual = (state, amplitude) -> [
        state[1]^2 - amplitude,
        state[2] - amplitude * state[1],
    ]
    energy = (state, amplitude) -> amplitude * state[1] + state[2]^2
    # Along the positive branch, E(A) = A^(3/2) + A^3.
    for amplitude in (1.0, 4.0, 9.0)
        state = [sqrt(amplitude), amplitude * sqrt(amplitude)]
        result = implicit_energy_derivative(residual, energy, state, amplitude)
        @test result.state_derivative ≈ [1 / (2sqrt(amplitude)), 3sqrt(amplitude) / 2]
        @test result.energy_derivative ≈ 3sqrt(amplitude) / 2 + 3amplitude^2
    end
end

@testset "Multiple observables include all amplitude dependence" begin
    residual = (state, amplitude) -> [
        state[1]^2 - amplitude,
        state[2] - amplitude * state[1],
    ]

    # Along this branch: u₁ = sqrt(A), u₂ = A^(3/2).
    observables = (state, amplitude) -> [
        amplitude * state[1] + state[2]^2, # A^(3/2) + A^3
        -state[1] * state[2],             # -A^2
        amplitude^2,                      # A^2
        state[2] - amplitude * state[1],  # zero
    ]

    for amplitude in (1.0, 4.0, 9.0)
        state = [sqrt(amplitude), amplitude * sqrt(amplitude)]
        result = implicit_observable_derivatives(
            residual, observables, state, amplitude,
        )

        @test result.state_derivative ≈ [
            1 / (2sqrt(amplitude)),
            3sqrt(amplitude) / 2,
        ]

        expected = [
            3sqrt(amplitude) / 2 + 3amplitude^2,
            -2amplitude,
            2amplitude,
            0.0,
        ]

        @test length(result.observable_derivatives) == length(expected)
        for (computed, reference) in zip(result.observable_derivatives, expected)
            @test computed ≈ reference atol=1e-12
        end
    end
end
@testset "Two-amplitude response and explicit observable dependence" begin
    # Exact positive branch: u1=sqrt(A1), u2=A2*u1.
    residual = (state, amplitudes) -> [state[1]^2-amplitudes[1],
        state[2]-amplitudes[2]*state[1]]
    observables = (state, amplitudes) -> [
        amplitudes[1]*state[1] + state[2]^2,
        state[1]*state[2] + amplitudes[2]^2]
    for amplitudes in ([4.0, 3.0], [0.01, 2.0])
        first, second = amplitudes
        state = [sqrt(first), second*sqrt(first)]
        expected_state = [1/(2sqrt(first)) 0; second/(2sqrt(first)) sqrt(first)]
        expected_observables = [1.5sqrt(first)+second^2 2first*second;
            second first+2second]
        result = implicit_observable_derivatives(residual, observables, state, amplitudes)
        @test result.state_derivative ≈ expected_state
        @test result.observable_derivatives ≈ expected_observables
        @test result.response_residual < 1e-14
        @test isfinite(result.state_jacobian_condition)
        @test implicit_branch_derivative(residual, state, amplitudes) ≈ expected_state
        @test state == [sqrt(first), second*sqrt(first)]
    end
    @test_throws DimensionMismatch implicit_branch_response(
        (state, amplitudes) -> [state[1]-amplitudes[1]], [1.0, 2.0], [1.0, 2.0])
    @test_throws SingularException implicit_branch_response(
        (state, amplitudes) -> [state[1]-amplitudes[1], 2state[1]-amplitudes[2]],
        [1.0, 0.0], [1.0, 2.0])
end
