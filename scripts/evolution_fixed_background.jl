# Run from the package root:  julia --project=. scripts/evolution_fixed_background.jl   (≈ 1 h)
# Long fixed-background evolution into the nonlinear regime (a = 0.9, α = 0.3, cosine).
using AxionCloud, Printf, DelimitedFiles
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
amplitudes = [1e-6, 1e-7]
span = (0.0, 1.5e9)
t = @elapsed result = evolve_fixed_background(problem, amplitudes, span; saveat=2.5e7, verbose=true)
println("done: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s"); flush(stdout)
rows = [(d.time, d.amplitudes..., real.(d.mode_frequencies)..., imag.(d.mode_frequencies)..., d.charges.energy, d.charges.angular_momentum, d.losses.energy, d.losses.angular_momentum)
        for d in result.diagnostics]
open(joinpath(@__DIR__, "..", "results", "fixed_background_a0.9_alpha0.3_diagnostics.csv"), "w") do io
    println(io, "time,A1,A2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,E,J,L_E,L_J")
    writedlm(io, [collect(row) for row in rows] |> x -> reduce(vcat, transpose.(x)), ',')
end
open(joinpath(@__DIR__, "..", "results", "fixed_background_a0.9_alpha0.3_solution.csv"), "w") do io
    println(io, "time,A1,A2")
    for (time, amps) in zip(result.ode_solution.t, result.ode_solution.u)
        @printf(io, "%.8e,%.8e,%.8e\n", time, amps[1], amps[2])
    end
end
flush(stdout)
