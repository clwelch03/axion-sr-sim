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
