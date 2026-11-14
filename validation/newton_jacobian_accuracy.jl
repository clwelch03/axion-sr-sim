# How accurate must the factorisation carried by a RadialSolution be for the evolution rates?
# At a converged state u* (problem p*), replace the factorisation by the one taken at the converged solution of a
# nearby problem p_ε (amplitudes × (1 + ε), b × (1 − 0.03ε), j × (1 − 0.1ε), roughly the direction of an evolution
# step), and compare the rates from balance_terms with those from the exact factorisation. "dist" is the max-norm
# distance between the two points (the measure solve_radial uses for `factorization_distance`), unscaled and in
# the scaled norm of the convergence test.
# Noise reference: the same state re-converged from p_ε's solution (unknowns differ at the Newton tolerance), with
# an exact factorisation at the re-converged point.
# States as in validation/newton_counts.jl (plateau, transition, nonlinear), built by amplitude continuation.
# Run from the package root (about 15 min; Slurm):  julia --project=. validation/newton_jacobian_accuracy.jl [states]
using AxionCloud, LinearAlgebra, Printf
const AC = AxionCloud

const DECAY_CONSTANT_RATIO = 0.0466873
const MODE_SET = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
const REFERENCE = let background = KerrBackground(0.9)
    RadialProblem(background, CosinePotential(), 0.3, MODE_SET; mesh=RadialMesh(background, 0.3, MODE_SET))
end
const STATES = Dict(
    "plateau" => ([0.00280773553, 3.3961332e-9], 0.9814753887857007, 0.8374944341429903),
    "transition" => ([0.0024024134667967956, 2.4071466062915234e-6], 0.9854426851242335, 0.8458838283206854),
    "nonlinear" => ([1.4719255751539253e-4, 4.100281105490423e-6], 0.9954602054965745, 0.8671482601144008),
)
problem_at(state::BackgroundState) =
    with_background(REFERENCE, KerrBackground(state.angular_momentum / state.mass^2), REFERENCE.mass_coupling * state.mass)

"""The RadialSolution at given unknowns with a fresh factorisation exactly there."""
function solution_from_unknowns(problem, amplitudes, unknowns, scale)
    layout = AC.UnknownLayout(problem)
    frequencies = AC._complexify(unknowns[AC.parameter_range(layout)][AC.frequency_range(layout)])
    scaling = AC.NewtonScaling(problem, layout, frequencies)
    factorization = AC.ScaledFactorization(AC.jacobian(problem, layout, unknowns, amplitudes, scale), scaling)
    residual_norm = norm(scaling.row .* AC.residual(problem, layout, unknowns, amplitudes, scale), Inf)
    parameters = AC.RadialParameters(layout, unknowns[AC.parameter_range(layout)], amplitudes, scale)
    RadialSolution(problem, layout, collect(amplitudes), float(scale), collect(unknowns), parameters, factorization, 0,
                   residual_norm, true)
end

"""Converged state by amplitude continuation on its background, polished at tol 1e-12 where reachable."""
function state_solution(amplitudes, background)
    problem = problem_at(background)
    solution = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=30))
        solution = solve_radial(problem, amplitudes .* factor; guess=solution)
        solution.converged || error("continuation failed at factor $factor")
    end
    polished = solve_radial(problem, amplitudes; guess=solution, tol=1e-12)
    best = polished.converged ? polished : solution
    return solution_from_unknowns(problem, amplitudes, best.unknowns, best.scale)
end

state_names = length(ARGS) ≥ 1 ? split(ARGS[1], ',') : ["plateau", "transition", "nonlinear"]
epsilons = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

function rates_with(solution, background, factorization)
    replaced = RadialSolution(solution.problem, solution.layout, solution.amplitudes, solution.scale, solution.unknowns,
                              solution.parameters, factorization, solution.iterations, solution.residual_norm, true)
    cache = EvolutionCache(REFERENCE; decay_constant_ratio=DECAY_CONSTANT_RATIO)
    terms = balance_terms(cache, replaced, background)
    return vcat(amplitude_rates(terms, cache) ./ solution.amplitudes, terms.state_rate.mass, terms.state_rate.angular_momentum)
end

for name in state_names
    amplitudes, mass, angular_momentum = STATES[name]
    background = BackgroundState(mass, angular_momentum)
    solution = state_solution(amplitudes, background)
    exact = rates_with(solution, background, solution.jacobian_factorization)
    println("\n=== $name: exact rates [dlnA1, dlnA2, db, dj]/dτ = ", exact)
    column = solution.jacobian_factorization.scaling.column
    rescaled_distance(other) = (u = AC.warm_start_unknowns(other, solution.scale); u .- solution.unknowns)
    @printf("%8s %11s %11s | %11s %11s %11s %11s\n", "eps", "dist(max)", "dist(scaled)", "err dlnA1", "err dlnA2", "err db", "err dj")
    for epsilon in epsilons
        perturbed_background = BackgroundState(mass * (1 - 0.03epsilon), angular_momentum * (1 - 0.1epsilon))
        nearby_problem, nearby_amplitudes = problem_at(perturbed_background), solution.amplitudes .* (1 + epsilon)
        nearby = solve_radial(nearby_problem, nearby_amplitudes; guess=solution, tol=1e-12)
        # tol 1e-12 is below the residual floor of the nonlinear state
        nearby.converged || (nearby = solve_radial(nearby_problem, nearby_amplitudes; guess=solution))
        nearby.converged || (println("eps $epsilon: nearby solve failed"); continue)
        nearby_exact = solution_from_unknowns(nearby.problem, nearby.amplitudes, nearby.unknowns, nearby.scale)
        difference = rescaled_distance(nearby_exact)
        errors = rates_with(solution, background, nearby_exact.jacobian_factorization) .- exact
        @printf("%8.0e %11.3e %11.3e | %11.3e %11.3e %11.3e %11.3e\n", epsilon, norm(difference, Inf),
                norm(difference ./ column, Inf), errors...)
        if epsilon == 1e-6
            # noise reference: re-converge at p* from p_ε's solution, exact factorisation there
            back = solve_radial(solution.problem, solution.amplitudes; guess=nearby_exact, tol=1e-10)
            back_exact = solution_from_unknowns(back.problem, back.amplitudes, back.unknowns, back.scale)
            noise = rates_with(back_exact, background, back_exact.jacobian_factorization) .- exact
            @printf("%8s %11.3e %11.3e | %11.3e %11.3e %11.3e %11.3e   (noise: re-converged at tol 1e-10)\n", "noise",
                    norm(back.unknowns .- solution.unknowns, Inf), norm((back.unknowns .- solution.unknowns) ./ column, Inf), noise...)
        end
        flush(stdout)
    end
end
