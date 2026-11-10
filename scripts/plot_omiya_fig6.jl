# Our joint evolutions on the axes of Omiya et al. 2211.01949 Fig. 6: energy per mode E_i/(f_a²M) against ω₁,I t,
# a = 0.99, α = 0.42 (left) and 0.3 (right), cosine (solid) and chiral (dashed) potentials, from the diagnostics CSVs of
# scripts/evolution_omiya_fig6.jl. Energy per mode is ω_i N_i (the energy carried by the quanta of phase i; their E_i
# is each cloud's own energy, Eq. 64). Missing runs are skipped. Repo palette (main: src/common/plotstyle.py).
# Dotted: Omiya's non-perturbative method (Sec. III: nonlinear primary, secondary on its averaged background, no
# deformation of the primary by the secondary), as reproduced on `main` from its coefficient tables
# (runs/adiabatic/multimode_secondary/results/two_cloud_evolution_alpha0{42,30}.csv: time in units of ω₁,I, E_i in
# f_a²M). That reproduction matched the paper's Fig. 6 at reading precision at α = 0.42; at α = 0.3 its peak was ≈ 10 %
# below and its plateaus 5–10 % above the paper's plotted values (that campaign's README).
# Run from the package root:  julia --project=scripts scripts/plot_omiya_fig6.jl [output.png] [omiya_results_dir]
using CairoMakie, DelimitedFiles, Printf

output_path = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "plots", "omiya_fig6_comparison.png")
results = joinpath(@__DIR__, "..", "results")
omiya_results = length(ARGS) ≥ 2 ? ARGS[2] :
                "/global/u1/c/clwelch/axion-sr-sim/runs/adiabatic/multimode_secondary/results"

function read_omiya(alpha)
    path = joinpath(omiya_results, @sprintf("two_cloud_evolution_alpha%03d.csv", round(Int, 100alpha)))
    isfile(path) || return nothing
    table, _ = readdlm(path, ','; header=true)
    return (; phase_time=table[:, 1], energy_211=table[:, 2], energy_322=table[:, 3])
end

function read_run(alpha, potential)
    path = joinpath(results, @sprintf("omiya_fig6_a0.99_alpha%.2f_%s_diagnostics.csv", alpha, potential))
    isfile(path) || return nothing
    lines = readlines(path)
    rate = parse(Float64, split(lines[1], "=")[2])
    table, header = readdlm(IOBuffer(join(lines[2:end], "\n")), ','; header=true)
    column(name) = table[:, findfirst(==(name), vec(header))]
    order = sortperm(column("time"))
    times = column("time")[order]
    keep = [index == 1 || times[index] > times[index - 1] for index in eachindex(times)]
    pick(name) = column(name)[order][keep]
    return (; phase_time=rate .* times[keep], energy_211=pick("Re_omega1") .* pick("N1"),
            energy_322=pick("Re_omega2") .* pick("N2"), total=pick("E"))
end

"""Log-linear interpolation of `values(times)` at `query_times` (both time grids increasing; queries inside the range)."""
function interpolate_log(times, values, query_times)
    map(query_times) do query
        upper = clamp(searchsortedfirst(times, query), 2, length(times))
        weight = (query - times[upper - 1]) / (times[upper] - times[upper - 1])
        exp((1 - weight) * log(values[upper - 1]) + weight * log(values[upper]))
    end
end

mode_211_color, mode_322_color = colorant"#3DAE2B", colorant"#00A0DF"
figure = Figure(size=(1300, 760), fontsize=15)
grid = (; xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06))
for (column_index, alpha) in enumerate((0.42, 0.3))
    axis = Axis(figure[1, column_index]; yscale=log10, title=@sprintf("a = 0.99,  α = %.2f", alpha),
                ylabel=column_index == 1 ? L"E_i / (f_a^2 M)" : "", grid...)
    difference_axis = Axis(figure[2, column_index]; xlabel=L"\omega_{1,I}\, t",
                           ylabel=column_index == 1 ? "ours (cosine) / Omiya − 1  [%]" : "", grid...)
    linkxaxes!(axis, difference_axis)
    hidexdecorations!(axis; grid=false)
    omiya = read_omiya(alpha)
    if omiya !== nothing
        lines!(axis, omiya.phase_time, omiya.energy_211; color=mode_211_color, linewidth=2.5, linestyle=:dot,
               label="|211⟩, Omiya method (cosine)")
        lines!(axis, omiya.phase_time, omiya.energy_322; color=mode_322_color, linewidth=2.5, linestyle=:dot,
               label="|322⟩, Omiya method (cosine)")
        peak = argmax(omiya.energy_211)
        @printf("α = %.2f Omiya : peak E₁ = %.1f at ω₁,I t = %.3f; final: E = (%.1f, %.1f)\n", alpha,
                omiya.energy_211[peak], omiya.phase_time[peak], omiya.energy_211[end], omiya.energy_322[end])
    end
    for (potential, style) in (("cosine", :solid), ("chiral", :dash))
        run = read_run(alpha, potential)
        run === nothing && continue
        lines!(axis, run.phase_time, run.energy_211; color=mode_211_color, linewidth=2, linestyle=style,
               label="|211⟩, $potential")
        lines!(axis, run.phase_time, run.energy_322; color=mode_322_color, linewidth=2, linestyle=style,
               label="|322⟩, $potential")
        peak = argmax(run.energy_211)
        @printf("α = %.2f %-6s: peak ω₂₁₁N₂₁₁ = %.1f at ω₁,I t = %.3f; final (ω₁,I t = %.2f): ωN = (%.1f, %.1f), E = %.1f\n",
                alpha, potential, run.energy_211[peak], run.phase_time[peak], run.phase_time[end], run.energy_211[end],
                run.energy_322[end], run.total[end])
        if potential == "cosine" && omiya !== nothing
            # Relative difference on Omiya's time grid, within the time range both runs cover.
            inside = run.phase_time[1] .≤ omiya.phase_time .≤ run.phase_time[end]
            query = omiya.phase_time[inside]
            for (ours, theirs, color) in ((run.energy_211, omiya.energy_211, mode_211_color),
                                          (run.energy_322, omiya.energy_322, mode_322_color))
                ratio = interpolate_log(run.phase_time, ours, query) ./ theirs[inside]
                lines!(difference_axis, query, 100 .* (ratio .- 1); color, linewidth=2)
                @printf("    %s: ours/Omiya − 1 ∈ [%+.2f, %+.2f] %%\n", color == mode_211_color ? "|211⟩" : "|322⟩",
                        100 * (minimum(ratio) - 1), 100 * (maximum(ratio) - 1))
            end
        end
    end
    hlines!(difference_axis, 0; color=(:black, 0.4), linewidth=1)
    xlims!(axis, 0, 10)
    ylims!(axis, 1, 3000)
end
rowsize!(figure.layout, 2, Relative(0.3))
# One legend for both panels (the entries are the same), in the lower right of the left panel, which is empty.
isempty(content(figure[1, 1]).scene.plots) ||
    axislegend(content(figure[1, 1]); position=:rb, framevisible=true, backgroundcolor=(:white, 0.85), nbanks=2)
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path")
