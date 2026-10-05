# julia --project=. tests/adiabatic/check_two_cloud_splines.jl
# Step 6a checks (piece 6, docs/multimode-adiabatic-design.md).
#   1. The not-a-knot spline reproduces a cubic exactly (value and derivative),
#      on unevenly spaced nodes.
#   2. On a smooth non-polynomial function its error falls by 16 per halving of
#      the node spacing (fourth order).
#   3. The saved coefficient tables: splines in E₁ built from every other
#      point, evaluated at the omitted points, against the saved values; and
#      dω₂/dE₁ from that spline against the one from all points. (The error of
#      the all-point spline is smaller, by about 16 if it is fourth order.)
using Printf, Random
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/two_cloud_evolution.jl"))

println("1. Cubic reproduced exactly")
Random.seed!(2)
nodes = sort(rand(12)) .* 3 .- 1
cubic(x) = 1 - 2x + 0.5x^2 + 0.3x^3
cubic_derivative(x) = -2 + x + 0.9x^2
spline = cubic_spline(nodes, cubic.(nodes))
grid = range(first(nodes), last(nodes); length=1001)
@printf("   max |S - p| = %.1e, max |S' - p'| = %.1e\n",
        maximum(abs(evaluate_spline(spline, x).value - cubic(x)) for x in grid),
        maximum(abs(evaluate_spline(spline, x).derivative - cubic_derivative(x)) for x in grid))

println("2. Error on sin(x) over [0, 3] as nodes double")
errors = map((8, 16, 32, 64)) do count
    nodes = collect(range(0, 3; length=count))
    spline = cubic_spline(nodes, sin.(nodes))
    grid = range(0, 3; length=2001)
    (maximum(abs(evaluate_spline(spline, x).value - sin(x)) for x in grid),
     maximum(abs(evaluate_spline(spline, x).derivative - cos(x)) for x in grid))
end
for i in eachindex(errors)
    ratio = i > 1 ? @sprintf("  ratios %.1f, %.1f", errors[i-1][1] / errors[i][1],
                             errors[i-1][2] / errors[i][2]) : ""
    @printf("   %2d nodes: value %.2e, derivative %.2e%s\n", (8, 16, 32, 64)[i], errors[i]..., ratio)
end

println("3. Coefficient splines in E1: every other point vs the omitted points")
coefficient_names = (:primary_omega, :secondary_omega, :primary_flux_per_energy,
                     :secondary_flux_per_energy, :m0_flux_coefficient, :m3_flux_coefficient)
for (alpha, energy_cap) in ((0.42, 1000.0), (0.3, Inf))
    table = read_two_cloud_table(alpha)
    energy = table.primary_energy
    kept = collect(1:2:length(energy))
    last(kept) == length(energy) || push!(kept, length(energy))
    omitted = setdiff(1:length(energy), kept)
    in_range = [i for i in omitted if table.primary_energy[i] <= energy_cap]
    @printf("   Mμ = %.2f: %d points, %d omitted with E1 <= %g\n",
            alpha, length(energy), length(in_range), energy_cap)
    for name in coefficient_names
        values = getproperty(table, name)
        half = cubic_spline(energy[kept], values[kept])
        scale = maximum(abs, values[[i for i in eachindex(values) if table.primary_energy[i] <= energy_cap]])
        errors_in_range = [abs(evaluate_spline(half, energy[i]).value - values[i]) / scale
                           for i in in_range]
        errors_all = [abs(evaluate_spline(half, energy[i]).value - values[i]) / maximum(abs, values)
                      for i in omitted]
        @printf("     %-26s max error / max |value|: %.1e in range, %.1e over the branch\n",
                name, maximum(errors_in_range), maximum(errors_all))
    end
    # dω₂/dE₁ from the full and the half spline, in range.
    full = cubic_spline(energy, table.secondary_omega)
    half = cubic_spline(energy[kept], table.secondary_omega[kept])
    slope(spline, i) = evaluate_spline(spline, energy[i]).derivative
    slopes = [slope(full, i) for i in in_range]
    differences = [abs(slope(half, i) - slope(full, i)) for i in in_range]
    @printf("     dω2/dE1: range %.2e .. %.2e, max |half - full| / max |full| = %.1e\n",
            minimum(slopes), maximum(slopes), maximum(differences) / maximum(abs, slopes))
end
