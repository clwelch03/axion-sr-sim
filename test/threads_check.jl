# Helper of the "Threaded kernels" test set in runtests.jl, run as a subprocess with a given thread count:
#   julia --project=<test environment> -t N threads_check.jl <output path>
# Writes, as raw Float64 values: a flag (1.0 if concurrent spindown evaluations on one cache and one guess equal the
# sequential ones bitwise, else 0.0), then the residual, the Jacobian values, the rates of every evaluation and the
# integrated densities with dual unknowns. The test compares the files of the 1-thread and 4-thread runs bitwise.
using AxionCloud, ForwardDiff
const AC = AxionCloud

background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
solution = solve_radial(problem, [3e-3, 1e-4])
layout, unknowns, amplitudes, scale = solution.layout, solution.unknowns, solution.amplitudes, solution.scale
cache = EvolutionCache(problem; decay_constant_ratio=1e-2, record_metric_work=true)
state = BackgroundState(1.0, 0.9)
terms = balance_terms(cache, solution, state)
tangents = SolutionTangents(solution, state, terms)
base_state = [log(3e-3), log(1e-4), 1.0, 0.9]
perturbations = [[3e-4, 0, 0, 0], [0, 3e-4, 0, 0], [0, 0, 1e-6, 0], [0, 0, 0, 1e-6], [-3e-4, 0, 0, 0], [0, 0, 0, -1e-6]]
sequential = [AC._spindown_evaluation(base_state .+ perturbation, cache, solution, tangents) for perturbation in perturbations]
concurrent = fetch.([Threads.@spawn AC._spindown_evaluation(base_state .+ perturbation, cache, solution, tangents)
                     for perturbation in perturbations])
same = all(index -> sequential[index].rates == concurrent[index].rates &&
                    sequential[index].solution.unknowns == concurrent[index].solution.unknowns &&
                    sequential[index].terms.action_jacobian == concurrent[index].terms.action_jacobian, eachindex(sequential))
dual_unknowns = [ForwardDiff.Dual(value, 0.1sin(index), cos(index)) for (index, value) in enumerate(unknowns)]
densities = AC.integrated_densities(problem, layout, dual_unknowns, scale)
values = vcat(same ? 1.0 : 0.0,
              AC.residual(problem, layout, unknowns .* (1 .+ 1e-3 .* sin.(1:length(unknowns))), amplitudes, scale),
              AC.jacobian(problem, layout, unknowns, amplitudes, scale).nzval,
              reduce(vcat, [evaluation.rates for evaluation in sequential]),
              ForwardDiff.value.(densities), reduce(vcat, [collect(ForwardDiff.partials(density)) for density in densities]))
write(ARGS[1], values)
