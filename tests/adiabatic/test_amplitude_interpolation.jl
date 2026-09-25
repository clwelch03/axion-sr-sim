@testset "Hermite energy and derivative reproduce a cubic" begin
    left = (amplitude=1.0, energy=1.0, energy_derivative=3.0)
    right = (amplitude=2.0, energy=8.0, energy_derivative=12.0)
    for amplitude in (1.0, 1.25, 1.5, 1.75, 2.0)
        result = interpolate_branch_energy(amplitude, left, right)
        @test result.energy ≈ amplitude^3
        @test result.energy_derivative ≈ 3amplitude^2
    end
end

@testset "Tabulated amplitude rate and interval boundaries" begin
    left = (amplitude=1.0, energy=1.0, energy_derivative=2.0,
            horizon_flux=-3.0, outer_flux=1.0)
    right = (amplitude=2.0, energy=4.0, energy_derivative=4.0,
             horizon_flux=-3.0, outer_flux=1.0)
    branch = (; amplitudes=[1.0, 2.0], samples=[left, right])
    for amplitude in (1.0, 1.5, 2.0)
        @test amplitude_rate(amplitude, branch) ≈ 1 / amplitude
    end
    @test branch_interval(1.5, [1.0, 2.0, 4.0]) == (left_index=1, right_index=2)
    @test branch_interval(2.0, [1.0, 2.0, 4.0]) == (left_index=2, right_index=3)
    @test branch_interval(4.0, [1.0, 2.0, 4.0]) == (left_index=2, right_index=3)
    @test_throws DomainError amplitude_rate(0.9, branch)
    @test_throws DomainError amplitude_rate(2.1, branch)
    @test_throws ArgumentError branch_interval(1.0, [1.0, 1.0])
end
