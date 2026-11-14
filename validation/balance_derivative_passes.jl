# The single dual pass of `balance_terms` (explicit and implicit background derivatives and the amplitude Jacobians
# together) against the earlier central difference of the explicit partial and against four re-solved radial problems,
# and the action-balance rates with the energy pieces off and on. Linear point and a nonlinear point, a = 0.9, α = 0.3.
# Differences are max|Δ| / max|value| per field of `BackgroundDerivatives`.
# Run from the package root:  julia --project=. validation/balance_derivative_passes.jl   (a few minutes; Slurm)
using AxionCloud, Printf
const AC = AxionCloud

background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
state = BackgroundState(1.0, background.a)

fields(derivatives) = [collect(getfield(derivatives, field)) for field in fieldnames(BackgroundDerivatives)]
relative_difference(first, second) = [maximum(abs, a .- b) / maximum(abs, b) for (a, b) in zip(fields(first), fields(second))]

for amplitudes in ([1e-6, 1e-7], [1e-2, 1e-4])
    solution = solve_radial(problem, amplitudes)
    cache = EvolutionCache(problem; decay_constant_ratio=1e-2, record_metric_work=true, last_solution=solution)
    println("A = $amplitudes")
    exact = background_charge_derivatives(cache, solution, state)
    resolved = balance_terms(cache, solution, state; resolve_background=true).derivatives
    @printf("  exact vs re-solved (step 1e-4):  %s\n", join((@sprintf("%.1e", d) for d in relative_difference(exact, resolved)), " "))
    for relative_step in (1e-5, 1e-6, 1e-7)
        central = background_charge_derivatives(cache, solution, state; relative_step)
        @printf("  exact vs central step %.0e:       %s\n", relative_step, join((@sprintf("%.1e", d) for d in relative_difference(exact, central)), " "))
        @printf("  central step %.0e vs re-solved:   %s\n", relative_step, join((@sprintf("%.1e", d) for d in relative_difference(central, resolved)), " "))
    end
    println("  (fields: charges ∂/∂b, charges ∂/∂j, actions ∂/∂b, actions ∂/∂j, frequencies ∂/∂b, frequencies ∂/∂j)")

    action_cache = EvolutionCache(problem; decay_constant_ratio=1e-2, last_solution=solution)
    recorded_terms = balance_terms(cache, solution, state)
    default_terms = balance_terms(action_cache, solution, state)
    recorded_rates, default_rates = amplitude_rates(recorded_terms, cache), amplitude_rates(default_terms, action_cache)
    @printf("  action rates, energy pieces on vs off: identical = %s, max|Δ| = %.2e (|rates| = %.2e)\n",
            recorded_rates == default_rates, maximum(abs, recorded_rates .- default_rates), maximum(abs, default_rates))
    @printf("  background rates identical = %s; action Jacobian identical = %s; action derivatives identical = %s\n",
            recorded_terms.state_rate == default_terms.state_rate, recorded_terms.action_jacobian == default_terms.action_jacobian,
            recorded_terms.derivatives.actions_wrt_mass == default_terms.derivatives.actions_wrt_mass &&
            recorded_terms.derivatives.actions_wrt_angular_momentum == default_terms.derivatives.actions_wrt_angular_momentum)
    for (label, terms_cache) in (("energy pieces off", action_cache), ("energy pieces on ", cache))
        balance_terms(terms_cache, solution, state)
        seconds = minimum(@elapsed(balance_terms(terms_cache, solution, state)) for _ in 1:3)
        @printf("  balance_terms, %s: %.3f s\n", label, seconds)
    end
    flush(stdout)
end
