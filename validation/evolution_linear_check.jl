# Fixed-background evolution in the linear regime: dA/dt must equal γA and a short Tsit5 run must
# reproduce exponential growth. Run from the package root:  julia --project=. validation/evolution_linear_check.jl
using AxionCloud, Printf
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
amplitudes = [1e-6, 1e-7]
solution = solve_radial(problem, amplitudes)
rates = amplitude_rates(solution)
println("dA/dt = ", rates, "  γA = ", growth_rates(solution) .* amplitudes, "  relative mismatch = ", rates ./ (growth_rates(solution) .* amplitudes) .- 1)
horizon = 1 / growth_rates(solution)[1]
result = evolve_fixed_background(problem, amplitudes, (0.0, 3horizon); saveat=horizon)
for (time, amps) in zip(result.ode_solution.t, result.ode_solution.u)
    @printf("t = %.3e  A = %s   linear prediction = %s\n", time, amps, amplitudes .* exp.(growth_rates(solution) .* time))
end
