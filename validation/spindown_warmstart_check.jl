# Why did the a₀ = 0.9, f_a/M_Pl = 3e-3 spindown run (job 61095796) fail? Its last good radial solve was a trial stage
# at t = 6.4187e8 (state P) that the integrator then rejected; the retry at t = 6.2924e8 (state F) was warm-started
# from P's solution and did not converge. Q is the last accepted step (t = 6.2742e8). States copied from
# results/spindown_a0.9_alpha0.3_fa3e-3-61095796.out (lines 135, 136, 129).
# Tests, each a single solve_radial at the target with the given guess:
#   1. F from P (reproduces the failure?)
#   2. amplitude jump only: P background, F amplitudes, from P
#   3. background jump only: F background, P amplitudes, from P
#   4. F from P with a first-order predictor u + (du/dA) ΔA + (du/d(b,j)) Δ(b,j)
#   5. F from Q (warm start from the last accepted step instead of the last trial stage)
# Run from the package root:  julia --project=. validation/spindown_warmstart_check.jl   (a few minutes)
using AxionCloud, ForwardDiff, Printf

mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
initial_background = KerrBackground(0.9)
initial_mass_coupling = 0.3
template = RadialProblem(initial_background, CosinePotential(), initial_mass_coupling, mode_set;
                         mesh=RadialMesh(initial_background, initial_mass_coupling, mode_set))
problem_at(mass, angular_momentum) =
    with_background(template, KerrBackground(angular_momentum / mass^2), initial_mass_coupling * mass)

state_P = (amplitudes=[0.013627504495350289, 1.3502696219047228e-8], mass=0.998653327542109, angular_momentum=0.8954546192995005)
state_F = (amplitudes=[0.011137401861956385, 6.2059368302084035e-9], mass=0.9990656139214019, angular_momentum=0.8968462721795618)
state_Q = (amplitudes=[0.010824296016147988, 5.80264118087624e-9], mass=0.9991123033422119, angular_momentum=0.8970038653409119)

"""Converged solution at a state by amplitude continuation from the linear regime on that background."""
function solve_by_continuation(state; n_steps=30)
    problem = problem_at(state.mass, state.angular_momentum)
    solution = nothing
    for factor in exp.(range(log(1e-3), 0.0; length=n_steps))
        solution = solve_radial(problem, state.amplitudes .* factor; guess=solution)
        solution.converged || error("continuation failed at factor $factor")
    end
    return solution
end

report(label, solution) = @printf("%-44s converged = %-5s  iterations = %2d  |F| = %.2e  ω = %s\n", label,
                                  solution.converged, solution.iterations, solution.residual_norm,
                                  join((@sprintf("%.10f%+.3ei", real(w), imag(w)) for w in mode_frequencies(solution)), ", "))

t = @elapsed begin
    solution_P = solve_by_continuation(state_P)
    solution_Q = solve_by_continuation(state_Q)
end
@printf("reference solves by continuation: %.1f s\n", t)
report("P (by continuation)", solution_P)
report("Q (by continuation)", solution_Q)
flush(stdout)

problem_F = problem_at(state_F.mass, state_F.angular_momentum)
report("1. F from P", solve_radial(problem_F, state_F.amplitudes; guess=solution_P)); flush(stdout)
report("2. amplitude jump only (P background)",
       solve_radial(problem_at(state_P.mass, state_P.angular_momentum), state_F.amplitudes; guess=solution_P)); flush(stdout)
report("3. background jump only (P amplitudes)", solve_radial(problem_F, state_P.amplitudes; guess=solution_P)); flush(stdout)

# 4. first-order predictor from P, in P's unknown scaling.
amplitude_directions = amplitude_response(solution_P)
mass_change = state_F.mass - state_P.mass
angular_momentum_change = state_F.angular_momentum - state_P.angular_momentum
residual_along_background = ForwardDiff.derivative(0.0) do epsilon
    AxionCloud.residual(problem_at(state_P.mass + epsilon * mass_change, state_P.angular_momentum + epsilon * angular_momentum_change),
                        solution_P.layout, solution_P.unknowns, solution_P.amplitudes, solution_P.scale)
end
background_direction = -(solution_P.jacobian_factorization \ residual_along_background)
predicted_unknowns = solution_P.unknowns .+ amplitude_directions * (state_F.amplitudes .- state_P.amplitudes) .+ background_direction
predicted = RadialSolution(solution_P.problem, solution_P.layout, solution_P.amplitudes, solution_P.scale, predicted_unknowns,
                           solution_P.parameters, solution_P.jacobian_factorization, 0, NaN, false)
report("4. F from P with first-order predictor", solve_radial(problem_F, state_F.amplitudes; guess=predicted)); flush(stdout)
report("5. F from Q (last accepted step)", solve_radial(problem_F, state_F.amplitudes; guess=solution_Q)); flush(stdout)
