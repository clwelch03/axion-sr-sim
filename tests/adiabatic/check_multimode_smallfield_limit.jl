# julia --project=. tests/adiabatic/check_multimode_background.jl
# Check B: c0 = <1 - cos ϕ₁> on saved cosine branch points against its
# small-field approximation c0_quad = <ϕ₁²>/2. At small amplitude,
# c0/c0_quad - 1 ≈ -c0_quad/4 (m = 1 component only).
using Serialization, Printf
include(joinpath(@__DIR__, "../../src/adiabatic/load.jl"))

branch_path = joinpath(
    @__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls",
)
saved = deserialize(branch_path)

# Rebuild the angular tables used when the branch was solved.
maximum_l = maximum(first, saved.modes)
grid = angular_grid(saved.n_polar, saved.n_azimuth)
polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(maximum_l, polar, grid)

amplitudes = [entry.amplitude for entry in saved.continuation]
# The six smallest saved amplitudes, then the points nearest A1 = 1, 3, 6, 13.
selected_indices = vcat(
    1:6, [argmin(abs.(amplitudes .- target)) for target in (1.0, 3.0, 6.0, 13.0)],
)

@printf("%-11s %-12s %-11s %-11s %-14s %-12s\n",
        "A1", "res. change", "c0_quad", "c0", "c0/c0_quad-1", "predicted")

for index in selected_indices
    entry = saved.continuation[index]
    solution = entry.solution

    # Reintegrate the saved solution with its own radial settings.
    profiles = nonlinear_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        saved.spin, saved.alpha, saved.modes, couplings, polar, azimuthal, grid;
        potential=:cosine, save_profile=true, saved.options...,
    )
    # Should be at round-off: the same inputs reproduce the saved solve.
    residual_change =
        maximum(abs, profiles.residual - solution.residual) / entry.amplitude

    # Compare at the (radius, polar node) where c0_quad is largest.
    largest_quadratic = 0.0
    deficit_there = 0.0
    for branch in (profiles.from_horizon, profiles.from_outer),
        sample in axes(branch.states, 2)

        radial_values = Dict(
            mode => branch.states[k, sample] for (k, mode) in enumerate(saved.modes)
        )
        field = reconstruct_field(radial_values, polar, azimuthal, grid)
        quadratic = azimuthal_average(field .^ 2 ./ 2, grid)
        polar_index = argmax(quadratic)

        if quadratic[polar_index] > largest_quadratic
            largest_quadratic = quadratic[polar_index]
            deficit_there = averaged_mass_deficit(
                radial_values, polar, azimuthal, grid,
            )[polar_index]
        end
    end

    relative_difference = deficit_there / largest_quadratic - 1
    predicted_difference = -largest_quadratic / 4
    @printf("%-11.5g %-12.3g %-11.5g %-11.5g %-14.5g %-12.5g\n",
            entry.amplitude, residual_change, largest_quadratic,
            deficit_there, relative_difference, predicted_difference)
end
