# Step control of the evolution drivers in the linear regime (a = 0.9, α = 0.3, cosine, default mesh,
# A₀ = (1e-6, 1e-7)), where A(t) = A₀ e^{γt} is known:
#   1. noise of the amplitude rates from the Newton tolerance (cold vs warm-started solves);
#   2. three e-folds of |211⟩ with the old settings (state A, Tsit5's initial-step heuristic), with state A and the
#      physics-based first step, and with the new defaults (state ln A, physics-based first step).
# Run from the package root:  julia --project=. validation/evolution_step_control.jl
using AxionCloud, Printf
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
amplitudes = [1e-6, 1e-7]

cold = solve_radial(problem, amplitudes; tol=1e-10)
nearby = solve_radial(problem, amplitudes .* (1 + 1e-4); tol=1e-10)
warm = solve_radial(problem, amplitudes; guess=nearby, tol=1e-10)
logarithmic_rates(solution) = amplitude_rates(solution) ./ amplitudes
growth = growth_rates(cold)
@printf("Im ω = (%.6e, %.6e)\n", growth...)
@printf("d ln A/dt cold = (%.10e, %.10e)\n", logarithmic_rates(cold)...)
@printf("d ln A/dt warm = (%.10e, %.10e)\n", logarithmic_rates(warm)...)
@printf("noise |cold − warm| = (%.1e, %.1e); relative to Im ω₂₁₁: (%.1e, %.1e)\n",
        abs.(logarithmic_rates(cold) .- logarithmic_rates(warm))..., (abs.(logarithmic_rates(cold) .- logarithmic_rates(warm)) ./ growth[1])...)
flush(stdout)

final_time = 3 / growth[1]
exact(time) = amplitudes .* exp.(growth .* time)
for (label, options) in (("old: state A, Tsit5 heuristic first step", (log_amplitudes=false, initial_step_fraction=nothing)),
                         ("state A, physics-based first step", (log_amplitudes=false,)),
                         ("new default: state ln A, physics-based first step", NamedTuple()))
    elapsed = @elapsed result = evolve_fixed_background(problem, amplitudes, (0.0, final_time); saveat=final_time / 3, options...)
    trajectory = amplitude_trajectory(result)
    error_211 = maximum(abs(trajectory.amplitudes[i, 1] / exact(trajectory.times[i])[1] - 1) for i in eachindex(trajectory.times))
    error_322 = maximum(abs(trajectory.amplitudes[i, 2] / exact(trajectory.times[i])[2] - 1) for i in eachindex(trajectory.times))
    @printf("%-52s %4d rhs evaluations, %6.1f s, max relative error vs A₀e^{γt}: 211 %.1e, 322 %.1e\n",
            label, length(result.diagnostics), elapsed, error_211, error_322)
    flush(stdout)
end
