# Plot a fixed-background evolution from the diagnostics CSV written by an evolution script:
#   julia --project=scripts scripts/plot_evolution.jl path/to/diagnostics.csv [output.png]
# Columns: time, A1, A2, Re_omega1, Re_omega2, Im_omega1, Im_omega2, E, J, L_E, L_J.
using CairoMakie, DelimitedFiles

diagnostics_path = length(ARGS) ≥ 1 ? ARGS[1] : error("give the diagnostics CSV path")
output_path = length(ARGS) ≥ 2 ? ARGS[2] : replace(diagnostics_path, ".csv" => ".png")
table, header = readdlm(diagnostics_path, ','; header=true)
column(name) = table[:, findfirst(==(name), vec(header))]
order = sortperm(column("time"))
time = column("time")[order]
positive = time .> 0

figure = Figure(size=(1300, 900))
amplitude_axis = Axis(figure[1, 1]; xlabel="t/M", ylabel="A", yscale=log10, title="Amplitudes")
lines!(amplitude_axis, time[positive], column("A1")[order][positive]; label="A₂₁₁")
lines!(amplitude_axis, time[positive], column("A2")[order][positive]; label="A₃₂₂")
axislegend(amplitude_axis; position=:lt)

energy_axis = Axis(figure[1, 2]; xlabel="t/M", ylabel="E / f_a²M", yscale=log10, title="Cloud energy and losses")
lines!(energy_axis, time[positive], column("E")[order][positive]; label="E")
losses = column("L_E")[order][positive]
lines!(energy_axis, time[positive][losses .> 0], losses[losses .> 0]; label="L_E > 0 (loss)", linestyle=:dash)
lines!(energy_axis, time[positive][losses .< 0], -losses[losses .< 0]; label="−L_E (gain)", linestyle=:dot)
axislegend(energy_axis; position=:lt)

rate_axis = Axis(figure[2, 1]; xlabel="t/M", ylabel="Im ω", title="Mode growth rates (eigenvalues)")
lines!(rate_axis, time[positive], column("Im_omega1")[order][positive]; label="Im ω₂₁₁")
lines!(rate_axis, time[positive], column("Im_omega2")[order][positive]; label="Im ω₃₂₂")
axislegend(rate_axis; position=:lt)

frequency_axis = Axis(figure[2, 2]; xlabel="t/M", ylabel="Re ω − Re ω(0)", title="Frequency shifts")
lines!(frequency_axis, time[positive], (column("Re_omega1")[order] .- column("Re_omega1")[order][1])[positive]; label="ω₂₁₁")
lines!(frequency_axis, time[positive], (column("Re_omega2")[order] .- column("Re_omega2")[order][1])[positive]; label="ω₃₂₂")
axislegend(frequency_axis; position=:lb)
save(output_path, figure)
println("saved ", output_path)
