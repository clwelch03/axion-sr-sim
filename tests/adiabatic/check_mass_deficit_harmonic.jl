# julia --project=. tests/adiabatic/check_mass_deficit_harmonic.jl
# Checks of mass_deficit_harmonic (piece 4a, docs/multimode-piece4-guide.md) on
# the saved Mμ = 0.42 cosine primary:
#   1. azimuthal_number = 0 reproduces averaged_mass_deficit;
#   2. the k = 2 harmonic tends to the small-field source (1/2) P1², with
#      P1 = Σ_l R_l1 Y_l1, and relative difference ∝ A1² (constant last column);
#   3. doubling the azimuthal grid changes the k = 2 harmonic only at round-off.
using Printf
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/continue_secondary.jl"))

primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
setup = load_secondary_setup(primary_path)
(; primary, grid, polar, azimuthal) = setup
match = primary.options.match

grid_fine = angular_grid(primary.n_polar, 2primary.n_azimuth)
azimuthal_fine = build_azimuthal_values(maximum(keys(azimuthal)), grid_fine.azimuth_nodes)

# The primary's R_lm at a tortoise radius, from the saved branch covering it.
function primary_values(profiles, position)
    branch = position <= match ? profiles.from_horizon : profiles.from_outer
    return primary_values_at(branch, position, primary.modes)
end

println("Checks 1 and 3 (maxima over polar nodes)")
for index in (1, 100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    for position in (-20.0, 5.0, 20.0, 50.0)
        values = primary_values(profiles, position)
        k0 = mass_deficit_harmonic(values, polar, azimuthal, grid, 0)
        average = averaged_mass_deficit(values, polar, azimuthal, grid)
        k2 = mass_deficit_harmonic(values, polar, azimuthal, grid, 2)
        k2_fine = mass_deficit_harmonic(values, polar, azimuthal_fine, grid_fine, 2)
        @printf("A1=%-8.4g r*=%5.1f  |k0 - average|=%.1e  |k2(2n) - k2(n)|/|k2|=%.1e\n",
                primary.continuation[index].amplitude, position,
                maximum(abs, k0 - average), maximum(abs, k2_fine - k2) / maximum(abs, k2))
    end
end

println("\nCheck 2: k = 2 harmonic against (1/2) P1² at r* = 20")
for index in 1:6
    (; profiles) = primary_profiles_at(setup, index)
    values = primary_values(profiles, 20.0)
    amplitude = primary.continuation[index].amplitude
    k2 = mass_deficit_harmonic(values, polar, azimuthal, grid, 2)
    small_field = sum(value .* polar[mode] for (mode, value) in values if mode[2] == 1) .^ 2 ./ 2
    relative = maximum(abs, k2 - small_field) / maximum(abs, small_field)
    @printf("A1=%-10.5g  relative difference=%.4e  /A1²=%.4f\n",
            amplitude, relative, relative / amplitude^2)
end
