# Parameter-space maps over (α, a) from scripts/plateau_scan.jl (results/plateau_scan_alpha*.csv), cosine potential,
# fixed background, |211⟩ + |322⟩ quasi-equilibrium (plateau). Panels:
#   1. plateau energy in |322⟩, ω₃₂₂ N₃₂₂ / (f_a²M);
#   2. rate-mismatch measure at the plateau, max_i |L_N,i^Omiya| / G_i (G_i: gross action flux through phase i);
#   3. first-order shift of Omiya's plateau relative to ours, max_i |δN_i / N_i|.
# Grey: |211⟩ not superradiant (no two-cloud plateau sought). ×: plateau Newton or the Omiya-type matching failed.
#
# Run from the package root:  julia --project=scripts scripts/plot_plateau_scan.jl [output.png]
using CairoMakie, DelimitedFiles, Printf

output_path = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "plots", "plateau_scan.png")
results = joinpath(@__DIR__, "..", "results")
paths = filter(name -> startswith(name, "plateau_scan_alpha") && endswith(name, ".csv"), readdir(results))
rows = Vector{Float64}[]
header = nothing
for path in paths
    table, file_header = readdlm(joinpath(results, path), ','; header=true)
    global header = vec(file_header)
    append!(rows, [Float64.(table[index, :]) for index in 1:size(table, 1)])
end
table = reduce(vcat, permutedims.(rows))
column(name) = table[:, findfirst(==(name), header)]
alphas, spins, status = column("alpha"), column("spin"), column("converged")
alpha_values, spin_values = sort(unique(alphas)), sort(unique(spins))

panels = [("plateau ω₃₂₂N₃₂₂ / (f_a²M)", column("E322"), (0, 4), v -> @sprintf("%.0f", v)),
          ("rate mismatch at plateau [%]", 100 .* column("measure"), (-2, 2), v -> @sprintf("%.2g", v)),
          ("Omiya plateau shift max|δN/N| [%]", 100 .* max.(abs.(column("shift_N1_rel")), abs.(column("shift_N2_rel"))),
           (-2, 2), v -> @sprintf("%.2g", v))]

mean(range) = (range[1] + range[2]) / 2
figure = Figure(size=(1500, 520), fontsize=14)
for (panel_index, (title, values, color_range, format)) in enumerate(panels)
    axis = Axis(figure[1, 2panel_index - 1]; title, xlabel="α = μM", ylabel=panel_index == 1 ? "a / M" : "")
    grid_values = fill(NaN, length(alpha_values), length(spin_values))
    for index in eachindex(values)
        i, j = findfirst(==(alphas[index]), alpha_values), findfirst(==(spins[index]), spin_values)
        status[index] == 1 && (grid_values[i, j] = log10(max(values[index], 1e-6)))
    end
    plot = heatmap!(axis, 1:length(alpha_values), 1:length(spin_values), grid_values; colormap=:viridis,
                    colorrange=color_range, nan_color=(:grey, 0.25))
    axis.xticks = (1:length(alpha_values), [@sprintf("%.2f", a) for a in alpha_values])
    axis.yticks = (1:length(spin_values), [@sprintf("%.2f", s) for s in spin_values])
    for index in eachindex(values)
        i, j = findfirst(==(alphas[index]), alpha_values), findfirst(==(spins[index]), spin_values)
        if status[index] == 1
            text!(axis, i, j; text=format(values[index]), align=(:center, :center), fontsize=11,
                  color=grid_values[i, j] > mean(color_range) ? :black : :white)
        elseif status[index] == 0
            text!(axis, i, j; text="×", align=(:center, :center), fontsize=16, color=:black)
        end
    end
    Colorbar(figure[1, 2panel_index], plot; label="log₁₀")
end
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path")
for index in sortperm(collect(zip(alphas, spins)))
    @printf("α = %.2f  a = %.2f  status %+d  plateau ωN = (%8.4g, %8.4g)  measure = %.3e  shift = (%+.3e, %+.3e)\n",
            alphas[index], spins[index], status[index], column("E211")[index], column("E322")[index],
            column("measure")[index], column("shift_N1_rel")[index], column("shift_N2_rel")[index])
end
