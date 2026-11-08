using AxionCloud, Printf, ForwardDiff
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
state = BackgroundState(1.0, background.a)
amplitudes = [1e-6, 1e-7]
solution = solve_radial(problem, amplitudes)
cache = AC.EvolutionCache(problem, problem.mass_coupling, 1e-4, false, solution, AC.EvolutionDiagnostics{Float64}[], 1e-10, false)
println("IFT step scan (∂E/∂b):")
for step in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
    wrt_mass, _ = background_charge_derivatives(cache, solution, state; relative_step=step)
    @printf("  step %.0e: %.10e\n", step, wrt_mass.energy)
end
println("re-solve FD step scan (∂E/∂b):")
for step in (1e-3, 1e-4, 1e-5)
    terms = balance_terms(cache, solution, state; resolve_background=true, relative_step=step)
    @printf("  step %.0e: %.10e\n", step, terms.charges_wrt_mass.energy)
end
println("α-only derivative at fixed background, dE/dα (current units):")
layout, unknowns, scale = solution.layout, solution.unknowns, solution.scale
delta = 1e-5 * 0.3
problem_plus, problem_minus = with_potential(RadialProblem(background, CosinePotential(), 0.3 + delta, mode_set), CosinePotential()), RadialProblem(background, CosinePotential(), 0.3 - delta, mode_set)
residual_partial = (AC.residual(problem_plus, layout, unknowns, amplitudes, scale) .- AC.residual(problem_minus, layout, unknowns, amplitudes, scale)) ./ (2delta)
charges_partial = (cloud_charges(problem_plus, layout, unknowns, scale).energy - cloud_charges(problem_minus, layout, unknowns, scale).energy) / (2delta)
direction = -(solution.jacobian_factorization \ residual_partial)
via_unknowns = ForwardDiff.derivative(e -> cloud_charges(problem, layout, unknowns .+ e .* direction, scale).energy, 0.0)
@printf("  IFT: explicit %.6e + response %.6e = %.10e\n", charges_partial, via_unknowns, charges_partial + via_unknowns)
E(α) = cloud_charges(solve_radial(RadialProblem(background, CosinePotential(), α, mode_set), amplitudes; guess=solution)).energy
@printf("  FD : %.10e\n", (E(0.3 + delta) - E(0.3 - delta)) / (2delta))
println("same for the spin direction at fixed α (should agree, as ∂Q/∂j did):")
spin_delta = 1e-5
bp, bm = KerrBackground(0.9 + spin_delta), KerrBackground(0.9 - spin_delta)
pp, pm = AC.with_background(problem, bp, 0.3), AC.with_background(problem, bm, 0.3)
residual_partial = (AC.residual(pp, layout, unknowns, amplitudes, scale) .- AC.residual(pm, layout, unknowns, amplitudes, scale)) ./ (2spin_delta)
charges_partial = (cloud_charges(pp, layout, unknowns, scale).energy - cloud_charges(pm, layout, unknowns, scale).energy) / (2spin_delta)
direction = -(solution.jacobian_factorization \ residual_partial)
via_unknowns = ForwardDiff.derivative(e -> cloud_charges(problem, layout, unknowns .+ e .* direction, scale).energy, 0.0)
@printf("  IFT: %.10e   FD: %.10e\n", charges_partial + via_unknowns,
        (cloud_charges(solve_radial(pp, amplitudes; guess=solution)).energy - cloud_charges(solve_radial(pm, amplitudes; guess=solution)).energy) / (2spin_delta))
