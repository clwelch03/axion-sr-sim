using Test
using Serialization

for filename in ("angular", "nonlinear", "tortoise", "radial",
                 "boundary_series", "boundaries", "integration",
                 "matching", "energy")
    include(joinpath(@__DIR__, "../../src/adiabatic", filename * ".jl"))
end

saved = deserialize(
    joinpath(@__DIR__, "../../runs/adiabatic/continuation.jls")
)
entry = first(saved.continuation)

grid = angular_grid(16, 31)
polar_values = build_polar_values(5, grid.cosine_theta_nodes)
azimuthal_values = build_azimuthal_values(5, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(5, polar_values, grid)

@testset "Energy: radial refinement at fixed endpoints" begin
    previous_solution = entry.solution
    energies = Float64[]

    for spacing in (0.05, 0.025, 0.0125)
        options = (
            step=spacing,
            inner=-100.0,
            outer=100.0,
            match=5.0,
        )

        println("Starting spacing = ", spacing)
        flush(stdout)

        # Re-establish matching at this spacing, using the previous
        # solution as the initial guess.
        solution = solve_nonlinear_mode(
            previous_solution.frequency, entry.amplitude,
            saved.spin, saved.alpha, saved.modes, couplings,
            polar_values, azimuthal_values, grid;
            horizon_guess=previous_solution.horizon_amplitudes,
            outer_tail_guess=previous_solution.outer_amplitudes[2:end],
            tolerance=1e-12,
            options...,
        )
        @test solution.scaled_residual_norm <= 1e-12

        profiles = nonlinear_matching_residual(
            solution.frequency,
            solution.horizon_amplitudes,
            solution.outer_amplitudes,
            saved.spin, saved.alpha, saved.modes, couplings,
            polar_values, azimuthal_values, grid;
            save_profile=true,
            options...,
        )

        energy = sum(
            integrate_energy_branch(
                branch, solution.frequency,
                saved.spin, saved.alpha, saved.modes,
                polar_values, azimuthal_values, grid,
            )
            for branch in (profiles.from_horizon, profiles.from_outer)
        )
        push!(energies, energy)

        println(
            "step=", spacing,
            ", E/Fa²=", energy,
            ", E/(Fa² A0²)=", energy / entry.amplitude^2,
            ", scaled residual=", solution.scaled_residual_norm,
        )
        flush(stdout)

        previous_solution = solution
    end

    coarse_difference = energies[1] - energies[2]
    fine_difference = energies[2] - energies[3]

    println("Coarse-minus-medium energy: ", coarse_difference)
    println("Medium-minus-fine energy: ", fine_difference)
    println(
        "Relative medium-to-fine change: ",
        abs(fine_difference / energies[3]),
    )

    if !iszero(fine_difference)
        println("Refinement ratio: ", coarse_difference / fine_difference)
    end
end