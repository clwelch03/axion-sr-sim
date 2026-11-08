using AxionCloud, Printf, LinearAlgebra
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
state = BackgroundState(1.0, background.a)
amplitudes = [1e-6, 1e-7]
solution = solve_radial(problem, amplitudes)
cache = AC.EvolutionCache(problem, problem.mass_coupling, 1e-4, false, solution, AC.EvolutionDiagnostics{Float64}[], 1e-10, false)
layout, unknowns, scale = solution.layout, solution.unknowns, solution.scale
function partial(step)
    plus, minus = BackgroundState(1 + step, state.angular_momentum), BackgroundState(1 - step, state.angular_momentum)
    (AC.residual(AC.current_problem(cache, plus), layout, unknowns, amplitudes, scale) .-
     AC.residual(AC.current_problem(cache, minus), layout, unknowns, amplitudes, scale)) ./ (2step)
end
reference = partial(1e-7)
n_block = AC.node_block_size(layout)
for step in (1e-4, 1e-5, 1e-6)
    difference = partial(step) .- reference
    worst = argmax(abs.(difference))
    row_block = (worst - 1) ÷ n_block + 1
    kind = row_block ≤ layout.n_nodes - 1 ? "interval $row_block (r* = $(problem.mesh.rstar[row_block]))" :
           row_block == layout.n_nodes ? "horizon rows" : "outer rows"
    @printf("step %.0e: max |Δ∂F/∂b| = %.3e at row %d (%s, local index %d); |∂F/∂b|_max = %.3e\n",
            step, abs(difference[worst]), worst, kind, (worst - 1) % n_block + 1, maximum(abs, reference))
    # how the direction du/db responds
    direction = solution.jacobian_factorization \ (partial(step) .- reference)
    worst_u = argmax(abs.(direction))
    u_block = (worst_u - 1) ÷ n_block + 1
    @printf("          max |Δ du/db| = %.3e at unknown %d (%s)\n", abs(direction[worst_u]), worst_u,
            u_block ≤ layout.n_nodes ? "node $u_block, r* = $(problem.mesh.rstar[u_block]), local $((worst_u - 1) % n_block + 1)" : "parameter block index $(worst_u - layout.n_nodes * n_block)")
end
# offsets and series at the first node vs b
for step in (0.0, 1e-5, 1e-4)
    plus = BackgroundState(1 + step, state.angular_momentum)
    p = AC.current_problem(cache, plus)
    @printf("b = 1 + %.0e: a = %.8f  s_1 = %.6e  s_end = %.6f  α = %.6f\n", step, p.background.a, p.mesh.offsets[1], p.mesh.offsets[end], p.mass_coupling)
end
