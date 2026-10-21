using Test
include(joinpath(@__DIR__, "../../../runs/adiabatic/multiphase_two_cloud/scripts/omiya_comparison.jl"))

@testset "Omiya comparison energy/amplitude conversion and physical time" begin
    growth1, growth2 = 0.03, 0.005
    initial_amplitudes = [1.0, 0.01]
    times = [0.0, 0.7, 1.9, 4.0]
    # The energies grow independently, while the unit-amplitude secondary
    # energy changes with the primary. A2 must include this normalization response.
    amplitudes = sort(unique(vcat(1.0, collect(range(0.7, 2.0; length=17)))))
    rows = [(; amplitude, primary_energy=4amplitude^2,
        secondary_energy=2 * (1 + 0.3 * 4amplitude^2),
        primary_omega=0.4, secondary_omega=0.415,
        primary_flux_per_energy=-2growth1, secondary_flux_per_energy=-2growth2,
        m0_flux_coefficient=0.0, m3_flux_coefficient=0.0) for amplitude in amplitudes]
    trajectory = omiya_trajectory(rows, initial_amplitudes, times)
    for (row, time) in zip(trajectory, times)
        expected_energy = 4exp(2growth1 * time)
        expected_secondary = initial_amplitudes[2] * exp(growth2 * time) *
            sqrt((1 + 0.3 * 4) / (1 + 0.3 * expected_energy))
        @test row.time == time
        @test row.primary_energy ≈ expected_energy rtol=1e-12
        @test row.A1 ≈ exp(growth1 * time) rtol=1e-12
        @test row.A2 ≈ expected_secondary rtol=1e-12
    end
    @test_throws ErrorException omiya_trajectory(rows, [0.91, 0.01], times)
    @test_throws DomainError omiya_trajectory(rows, initial_amplitudes, [0.0, 100.0])
end
