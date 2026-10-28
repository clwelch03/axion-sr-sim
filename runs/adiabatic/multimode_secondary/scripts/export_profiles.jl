# julia --project=. runs/adiabatic/multimode_secondary/scripts/export_profiles.jl [tag target_energy]
# Export equatorial <1 - cos ϕ₁> and the radial profiles |R_11| (primary) and
# |R_22| (averaged secondary) at a few primary amplitudes, for plot_secondary.py.
# Without arguments: the Mμ = 0.42 branch, target energy 1000, output work/profiles.csv.
# With a tag (e.g. alpha030): results/<tag>_cosine_branch.jls and
# results/<tag>_cosine_secondary.jls, output work/<tag>_profiles.csv (ignored).
using Serialization, Printf
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

campaign = joinpath(@__DIR__, "..")
tag = length(ARGS) >= 1 ? ARGS[1] : ""
target_energy = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1000.0
primary_path = isempty(tag) ?
    joinpath(@__DIR__, "../../fig2_20260922/results/cosine_branch.jls") :
    joinpath(campaign, "results/$(tag)_cosine_branch.jls")
secondary_path = joinpath(campaign, isempty(tag) ? "results/cosine_secondary.jls" :
                                                   "results/$(tag)_cosine_secondary.jls")
output_path = joinpath(campaign, isempty(tag) ? "work/profiles.csv" : "work/$(tag)_profiles.csv")
setup = load_primary_secondary_setup(primary_path)
secondary_branch = deserialize(secondary_path)
secondary_branch.reading == "averaged" || error("expected the averaged secondary branch")
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup

# A one-node angular "grid" at the equator, x = cos θ = 0, sharing the azimuthal nodes.
equator_grid = (; cosine_theta_nodes=[0.0], polar_weights=[1.0],
                  azimuth_nodes=grid.azimuth_nodes, azimuth_weight=grid.azimuth_weight)
equator_polar = build_polar_values(6, equator_grid.cosine_theta_nodes)

amplitudes = [entry.amplitude for entry in primary.continuation]
# Also the point whose primary energy is closest to target_energy (Fig. 7 analogue).
energies = [row.energy for row in primary.rows]
selected_indices = unique(vcat(
    [argmin(abs.(amplitudes .- target)) for target in (0.0, 0.1, 1.0, 3.0, 6.0, 13.0)],
    argmin(abs.(energies .- target_energy)),
))

mkpath(joinpath(campaign, "work"))
open(output_path, "w") do io
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
println("wrote ", output_path)
