# julia --project=. tests/adiabatic/primary_secondary/check_two_cloud_evolution.jl
# Step 6c checks of the Fig. 6 evolution (piece 6, docs/multimode-adiabatic-design.md),
# on the features read off the figure: peak of E1 and E2 (value, time) and the
# energies at ω_{1,I} t = 10.
#   1. Time step 2e-3, 1e-3 (default), 5e-4: changes should fall by about 16.
#   2. Coefficient splines built from every other saved point: the change
#      bounds the interpolation error of the all-point splines.
#   3. Initial straight line: ln E1 grows with slope 2 in ω_{1,I} t while E2
#      is negligible.
using Printf
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

feature_names = (:primary_peak, :primary_peak_time, :secondary_peak,
                 :secondary_peak_time, :primary_final, :secondary_final)
print_features(label, features) = @printf("  %-28s %s\n", label,
    join([@sprintf("%.6g", getproperty(features, name)) for name in feature_names], "  "))
relative_changes(a, b) = join([@sprintf("%.1e", abs(getproperty(a, name) / getproperty(b, name) - 1))
                               for name in feature_names], "  ")

for alpha in (0.42, 0.3)
    @printf("Mμ = %.2f; columns: %s\n", alpha, join(string.(feature_names), ", "))
    reference = evolution_features(fig6_evolution(alpha))
    print_features("default (step 1e-3)", reference)

    coarse = evolution_features(fig6_evolution(alpha; time_step=2e-3))
    fine = evolution_features(fig6_evolution(alpha; time_step=5e-4))
    println("  1. relative change, step 2e-3 vs 1e-3: ", relative_changes(coarse, reference))
    println("     relative change, step 1e-3 vs 5e-4: ", relative_changes(reference, fine))

    count = length(read_two_cloud_table(alpha).primary_energy)
    points = collect(1:2:count)
    last(points) == count || push!(points, count)
    half = evolution_features(fig6_evolution(alpha; points))
    println("  2. relative change, every other point:  ", relative_changes(half, reference))

    evolution = fig6_evolution(alpha)
    early = findall(<=(1.0), evolution.times)
    slope = (log(evolution.primary_energies[early[end]]) - log(evolution.primary_energies[1])) /
            (evolution.times[early[end]] - evolution.times[1])
    @printf("  3. d ln E1 / d(ω_{1,I} t) over [0, 1]: %.6f; E2 at t = 1: %.2e\n",
            slope, evolution.secondary_energies[early[end]])
end
