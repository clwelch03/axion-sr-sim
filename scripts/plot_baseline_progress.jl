# Occupations and energies along a finished fixed-background evolution, from the diagnostics CSV written by
# scripts/evolution_fixed_background.jl (one row per rhs evaluation). Top: occupations N_i. Bottom: total cloud energy E
# and ω_i N_i, the energy carried by each phase's quanta (E = Σ ω_i N_i − L̄; the difference is the averaged
# Lagrangian, which starts at quartic order). Left: the whole run; right: the transition.
# Run from the package root:  julia --project=scripts scripts/plot_baseline_progress.jl [diagnostics.csv] [output.png]
using CairoMakie, DelimitedFiles, Printf

csv_path = length(ARGS) ≥ 1 ? ARGS[1] :
           joinpath(@__DIR__, "..", "results", "fixed_background_a0.9_alpha0.3_contained_diagnostics.csv")
output_path = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(@__DIR__, "..", "plots", "baseline_progress.png")

table, header = readdlm(csv_path, ','; header=true)
column(name) = table[:, findfirst(==(name), vec(header))]
order = sortperm(column("time"))               # rhs evaluations include integrator stages; sort and drop repeats
times = column("time")[order]
keep = [index == 1 || times[index] > times[index - 1] for index in eachindex(times)]
times = times[keep]
series(name) = column(name)[order][keep]
occupation_211, occupation_322 = series("N1"), series("N2")
energy = series("E")
phase_energy_211 = series("Re_omega1") .* occupation_211
phase_energy_322 = series("Re_omega2") .* occupation_322

# repo palette (main: src/common/plotstyle.py): brand colours in cycle order, C0 green, C1 blue
mode_211_color, mode_322_color, ink = colorant"#3DAE2B", colorant"#00A0DF", colorant"#0b0b0b"
figure = Figure(size=(1300, 820), fontsize=15)
Label(figure[0, 1:2], "Baseline evolution: a = 0.9, α = 0.3, cosine potential, fixed background, equal initial occupations";
      fontsize=17, color=ink, halign=:left, tellwidth=false)
windows = ((times[1], times[end], "Whole run"), (6.3e8, 7.5e8, "Transition"))
for (column_index, (time_min, time_max, title)) in enumerate(windows)
    shown = time_min .≤ times .≤ time_max
    axis_n = Axis(figure[1, column_index]; yscale=log10, title, titlealign=:left,
                  ylabel=column_index == 1 ? L"occupation $\mathcal{N}_i$  [$f_a^2 M^2$]" : "",
                  xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06))
    lines!(axis_n, times[shown], occupation_211[shown]; color=mode_211_color, linewidth=2, label=L"|211\rangle")
    lines!(axis_n, times[shown], occupation_322[shown]; color=mode_322_color, linewidth=2, label=L"|322\rangle")
    axis_e = Axis(figure[2, column_index]; yscale=log10, xlabel=L"t/M", ylabel=column_index == 1 ? L"energy  [$f_a^2 M$]" : "",
                  xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06))
    lines!(axis_e, times[shown], energy[shown]; color=ink, linewidth=2, label=L"total cloud energy $E$")
    lines!(axis_e, times[shown], phase_energy_211[shown]; color=mode_211_color, linewidth=2, linestyle=:dash,
           label=L"\omega_{211}\,\mathcal{N}_{211}")
    lines!(axis_e, times[shown], phase_energy_322[shown]; color=mode_322_color, linewidth=2, linestyle=:dash,
           label=L"\omega_{322}\,\mathcal{N}_{322}")
    xlims!(axis_n, time_min, time_max)
    xlims!(axis_e, time_min, time_max)
    ticks = column_index == 1 ? [0.0, 5e8, 1e9] : collect(6.4e8:0.2e8:7.4e8)
    axis_n.xticks = ticks
    axis_e.xticks = ticks
    if column_index == 1
        axislegend(axis_n; position=:lt, framevisible=false)
        axislegend(axis_e; position=:lt, framevisible=false)
    end
end
colgap!(figure.layout, 40)
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path ($(length(times)) points)")
peak = argmax(occupation_211)
@printf("peak 𝒩₂₁₁ = %.4e at t = %.4e (E = %.4e); final: 𝒩 = (%.4e, %.4e), E = %.4e, 𝒩₃₂₂/𝒩₂₁₁ = %.3f\n",
        occupation_211[peak], times[peak], energy[peak], occupation_211[end], occupation_322[end], energy[end],
        occupation_322[end] / occupation_211[end])
