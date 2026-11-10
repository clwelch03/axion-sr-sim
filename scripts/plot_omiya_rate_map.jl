# Map of the rate difference between the joint (mutually deforming) solution and Omiya-type rates (Omiya et al.
# 2211.01949 Sec. III: |211⟩ independent of |322⟩, |322⟩ linear on it), at equal occupations, from the column CSVs of
# scripts/omiya_rate_map.jl, with our Fig. 6 cosine trajectories (scripts/evolution_omiya_fig6.jl) drawn on top.
#
# Axes: energy carried by each phase, ω_i N_i in f_a²M (the trajectory coordinates of plots/omiya_fig6_comparison.png).
# Colour: max over the two phases of |L_N,joint − L_N,Omiya| / G_i, where G_i = Σ_q |q_i| |loss_q| is the gross action
# flux through phase i (joint solution). The plain ratio L_joint/L_Omiya − 1 is not used because L_N₂₁₁ changes sign
# on the grid. The numbers printed on the markers are that quantity in percent.
#
# Run from the package root:  julia --project=scripts scripts/plot_omiya_rate_map.jl [output.png]
using CairoMakie, DelimitedFiles, Printf

output_path = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "plots", "omiya_rate_map.png")
results = joinpath(@__DIR__, "..", "results")

function read_map(alpha)
    paths = filter(name -> startswith(name, @sprintf("omiya_rate_map_a0.99_alpha%.2f_col", alpha)) && endswith(name, ".csv"),
                   readdir(results))
    rows = Vector{Float64}[]
    header = nothing
    for path in paths
        table, file_header = readdlm(joinpath(results, path), ','; header=true)
        header = vec(file_header)
        size(table, 1) == 0 && continue
        append!(rows, [Float64.(table[index, :]) for index in 1:size(table, 1)])
    end
    table = reduce(vcat, permutedims.(rows))
    column(name) = table[:, findfirst(==(name), header)]
    labels = [name[length("joint_loss_q") + 1:end] for name in header if startswith(name, "joint_loss_q")]
    charges = [parse.(Int, split(label, "_")) for label in labels]
    joint_components = reduce(hcat, column("joint_loss_q" * label) for label in labels)
    gross = [sum(abs(charges[k][phase]) * abs.(joint_components[:, k]) for k in eachindex(charges)) for phase in 1:2]
    difference = [abs.(column("joint_L_N$phase") .- column("omiya_L_N$phase")) ./ gross[phase] for phase in 1:2]
    return (; energy_211=column("Re_omega1") .* column("N1"), energy_322=column("Re_omega2") .* column("N2"),
            target_211=column("target_E211"), target_322=column("target_E322"),
            measure=max.(difference[1], difference[2]))
end

function read_trajectory(alpha)
    path = joinpath(results, @sprintf("omiya_fig6_a0.99_alpha%.2f_cosine_diagnostics.csv", alpha))
    lines = readlines(path)
    table, header = readdlm(IOBuffer(join(lines[2:end], "\n")), ','; header=true)
    column(name) = table[:, findfirst(==(name), vec(header))]
    order = sortperm(column("time"))
    return (; energy_211=(column("Re_omega1") .* column("N1"))[order], energy_322=(column("Re_omega2") .* column("N2"))[order])
end

colormap = :viridis
color_range = (-4, 0)   # log₁₀ of the measure
figure = Figure(size=(1300, 620), fontsize=15)
for (column_index, alpha) in enumerate((0.42, 0.3))
    axis = Axis(figure[1, column_index]; xscale=log10, yscale=log10, title=@sprintf("a = 0.99,  α = %.2f,  cosine", alpha),
                xlabel=L"\omega_{211} \mathcal{N}_{211} \; / \; (f_a^2 M)",
                ylabel=column_index == 1 ? L"\omega_{322} \mathcal{N}_{322} \; / \; (f_a^2 M)" : "",
                xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06))
    map_data = read_map(alpha)
    trajectory = read_trajectory(alpha)
    visible = trajectory.energy_322 .> 0.3
    lines!(axis, trajectory.energy_211[visible], trajectory.energy_322[visible]; color=:black, linewidth=1.5,
           label="our Fig. 6 trajectory")
    scatter!(axis, [trajectory.energy_211[end]], [trajectory.energy_322[end]]; color=:black, marker=:star5,
             markersize=16, label="end (ω₁,I t = 10)")
    scatter!(axis, map_data.energy_211, map_data.energy_322; color=log10.(max.(map_data.measure, 1e-6)), colormap,
             colorrange=color_range, marker=:rect, markersize=34, strokecolor=:white, strokewidth=1)
    for index in eachindex(map_data.measure)
        percent = 100map_data.measure[index]
        text!(axis, map_data.energy_211[index], map_data.energy_322[index];
              text=percent ≥ 10 ? @sprintf("%.0f", percent) : percent ≥ 0.1 ? @sprintf("%.1f", percent) : "<.1",
              align=(:center, :center), fontsize=10, color=map_data.measure[index] > 0.05 ? :black : :white)
        @printf("α = %.2f  E = (%6g, %6g)  ωN = (%9.4g, %9.4g)  measure = %.3e\n", alpha, map_data.target_211[index],
                map_data.target_322[index], map_data.energy_211[index], map_data.energy_322[index], map_data.measure[index])
    end
    xlims!(axis, 10, 3e4)
    ylims!(axis, 0.3, 3e4)
    axislegend(axis; position=:lt, framevisible=true, backgroundcolor=(:white, 0.85))
end
Colorbar(figure[1, 3]; colormap, limits=color_range, label="log₁₀ max_i |ΔL_N,i| / gross flux_i")
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path")
