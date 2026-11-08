# Implicit-function-theorem background derivatives of the charges versus the four-re-solve finite
# differences, at the a = 0.9 initial point and at a nonlinear point.
# Run from the package root:  julia --project=. validation/background_derivatives.jl
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
state = BackgroundState(1.0, background.a)
for amplitudes in ([1e-6, 1e-7], [1e-2, 1e-4])
    solution = solve_radial(problem, amplitudes)
    cache = EvolutionCache(problem; decay_constant_ratio=1e-2, last_solution=solution)
    t_ift = @elapsed derivatives = background_charge_derivatives(cache, solution, state)
    t_ift = @elapsed derivatives = background_charge_derivatives(cache, solution, state)
    ift_mass, ift_spin = derivatives.charges_wrt_mass, derivatives.charges_wrt_angular_momentum
    t_fd = @elapsed terms = balance_terms(cache, solution, state; resolve_background=true)
    println("A = $amplitudes")
    @printf("  IFT (%.1fs): ∂Q/∂b = (%.10e, %.10e)  ∂Q/∂j = (%.10e, %.10e)\n", t_ift, ift_mass.energy, ift_mass.angular_momentum, ift_spin.energy, ift_spin.angular_momentum)
    @printf("  FD  (%.1fs): ∂Q/∂b = (%.10e, %.10e)  ∂Q/∂j = (%.10e, %.10e)\n", t_fd, terms.derivatives.charges_wrt_mass.energy, terms.derivatives.charges_wrt_mass.angular_momentum,
            terms.derivatives.charges_wrt_angular_momentum.energy, terms.derivatives.charges_wrt_angular_momentum.angular_momentum)
    @printf("  relative differences: %.2e %.2e %.2e %.2e\n", ift_mass.energy / terms.derivatives.charges_wrt_mass.energy - 1,
            ift_mass.angular_momentum / terms.derivatives.charges_wrt_mass.angular_momentum - 1,
            ift_spin.energy / terms.derivatives.charges_wrt_angular_momentum.energy - 1,
            ift_spin.angular_momentum / terms.derivatives.charges_wrt_angular_momentum.angular_momentum - 1)
    @printf("  actions:  IFT ∂N/∂b = %s  FD %s\n", derivatives.actions_wrt_mass, terms.derivatives.actions_wrt_mass)
    @printf("            IFT ∂N/∂j = %s  FD %s\n", derivatives.actions_wrt_angular_momentum, terms.derivatives.actions_wrt_angular_momentum)
    t_fast = @elapsed fast = balance_terms(cache, solution, state)
    @printf("  balance_terms with IFT: %.1fs; A' = %s vs %s\n", t_fast, amplitude_rates(fast), amplitude_rates(terms))
    flush(stdout)
end
