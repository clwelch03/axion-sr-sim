# julia --project=. runs/adiabatic/multimode_secondary/scripts/export_profiles.jl
# Export equatorial <1 - cos ϕ₁> and the radial profiles |R_11| (primary) and
# |R_22| (averaged secondary) at a few primary amplitudes, for plot_secondary.py.
# Output: runs/adiabatic/multimode_secondary/work/profiles.csv (ignored).
using Serialization, Printf
include(joinpath(@__DIR__, "../../../../src/adiabatic/drivers/continue_secondary.jl"))

campaign = joinpath(@__DIR__, "..")
setup = load_secondary_setup(joinpath(
    @__DIR__, "../../fig2_20260922/results/cosine_branch.jls"))
secondary_branch = deserialize(joinpath(campaign, "results/cosine_secondary.jls"))
secondary_branch.reading == "averaged" || error("expected the averaged secondary branch")
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup

# A one-node angular "grid" at the equator, x = cos θ = 0, sharing the azimuthal nodes.
equator_grid = (; cosine_theta_nodes=[0.0], polar_weights=[1.0],
                  azimuth_nodes=grid.azimuth_nodes, azimuth_weight=grid.azimuth_weight)
equator_polar = build_polar_values(6, equator_grid.cosine_theta_nodes)

amplitudes = [entry.amplitude for entry in primary.continuation]
# Also the point whose primary energy is closest to 1000 Fa² M (Fig. 7 analogue).
energies = [row.energy for row in primary.rows]
selected_indices = unique(vcat(
    [argmin(abs.(amplitudes .- target)) for target in (0.0, 0.1, 1.0, 3.0, 6.0, 13.0)],
    argmin(abs.(energies .- 1000)),
))

mkpath(joinpath(campaign, "work"))
open(joinpath(campaign, "work/profiles.csv"), "w") do io
    println(io, "amplitude,tortoise_radius,mass_deficit_equator,primary_abs_R11,secondary_abs_R22")
    for index in selected_indices
        (; profiles) = primary_profiles_at(setup, index)
        solution = secondary_branch.rows[index].solution
        secondary = secondary_matching_residual(
            solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
            profiles, primary.modes, primary.spin, primary.alpha, secondary_modes,
            couplings, polar, azimuthal, grid; save_profile=true, reading=:averaged,
        )
        # Pair each secondary branch with the primary branch it read from.
        branch_pairs = ((secondary.from_horizon, profiles.from_horizon),
                        (secondary.from_outer, profiles.from_outer))
        for (secondary_part, primary_part) in branch_pairs
            for (sample, position) in enumerate(secondary_part.positions)
                primary_values = primary_values_at(primary_part, position, primary.modes)
                deficit = averaged_mass_deficit(
                    primary_values, equator_polar, azimuthal, equator_grid)[1]
                @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g\n",
                        amplitudes[index], position, deficit,
                        abs(primary_values[(1, 1)]), abs(secondary_part.states[1, sample]))
            end
        end
    end
end
println("wrote ", joinpath(campaign, "work/profiles.csv"))
