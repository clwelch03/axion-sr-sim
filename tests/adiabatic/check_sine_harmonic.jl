# julia --project=. tests/adiabatic/check_sine_harmonic.jl
# Checks of sine_harmonic (piece 5a, docs/multimode-adiabatic-design.md) on the
# saved Mμ = 0.42 cosine primary:
#   1. sin ϕ₁ has only odd harmonics of ψ₁: k = 0 and k = 2 vanish;
#   2. the k = 1 harmonic tends to P1 = Σ_l R_l1 Y_l1, with relative difference
#      ∝ A1² (constant last column);
#   3. doubling the azimuthal grid changes the k = 1 harmonic only at round-off.
using Printf
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/continue_secondary.jl"))

primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
setup = load_secondary_setup(primary_path)
(; primary, grid, polar, azimuthal) = setup
match = primary.options.match

grid_fine = angular_grid(primary.n_polar, 2primary.n_azimuth)
azimuthal_fine = build_azimuthal_values(maximum(keys(azimuthal)), grid_fine.azimuth_nodes)
@printf("angular grid: %d polar x %d azimuthal\n", primary.n_polar, primary.n_azimuth)

# The primary's R_lm at a tortoise radius, from the saved branch covering it.
function primary_values(profiles, position)
    branch = position <= match ? profiles.from_horizon : profiles.from_outer
    return primary_values_at(branch, position, primary.modes)
end
first_harmonic_profile(values) =
    sum(value .* polar[mode] for (mode, value) in values if mode[2] == 1)

println("Checks 1 and 3 (maxima over polar nodes, relative to max |k1|)")
for index in (1, 100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    for position in (-20.0, 5.0, 20.0, 50.0)
        values = primary_values(profiles, position)
        k1 = sine_harmonic(values, polar, azimuthal, grid, 1)
        scale = maximum(abs, k1)
        k0 = sine_harmonic(values, polar, azimuthal, grid, 0)
        k2 = sine_harmonic(values, polar, azimuthal, grid, 2)
        k1_fine = sine_harmonic(values, polar, azimuthal_fine, grid_fine, 1)
        @printf("A1=%-8.4g r*=%5.1f  max|k1|=%.2e  |k0|=%.1e |k2|=%.1e  |k1(2n)-k1(n)|=%.1e  max|P1|=%.2e\n",
                primary.continuation[index].amplitude, position, scale,
                maximum(abs, k0) / scale, maximum(abs, k2) / scale,
                maximum(abs, k1_fine - k1) / scale, maximum(abs, first_harmonic_profile(values)))
    end
end

println("\nCheck 2: k = 1 harmonic against P1 at r* = 20")
for index in (1:6..., 100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    values = primary_values(profiles, 20.0)
    amplitude = primary.continuation[index].amplitude
    k1 = sine_harmonic(values, polar, azimuthal, grid, 1)
    small_field = first_harmonic_profile(values)
    relative = maximum(abs, k1 - small_field) / maximum(abs, small_field)
    @printf("A1=%-10.5g  relative difference=%.4e  /A1²=%.4f\n",
            amplitude, relative, relative / amplitude^2)
end
