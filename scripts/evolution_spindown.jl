# Run from the package root:
#   julia --project=. scripts/evolution_spindown.jl <f_a/M_Pl> [initial spin, default 0.9] [final t/M₀] [restart CSV]
# Restart: the 4th argument is the *_solution.csv of an earlier run with the same f_a and initial spin; the run
# continues from its last row (t, A, b, j) to the final time, in the same reference units and on the same mesh. The
# radial solution at that state is rebuilt by amplitude continuation (a cold solve can fail in the nonlinear regime).
# Restarted outputs get the suffix _restart.
# (1–2 h; submit to Slurm). The baseline of scripts/evolution_fixed_background.jl (a = 0.9 unless given, α = 0.3, cosine, contained mesh, equal initial
# occupations, t ≤ 1.5e9) with the black hole absorbing the horizon fluxes: mass b = M/M₀ and angular momentum
# j = J/M₀² evolve, coupled through (f_a/M_Pl)² (non-reduced Planck mass). Action balance (default), so no
# metric-work term enters; w is still recorded. Reference units are the initial mass M₀.
# The mesh is built for the initial spin and reused, in units of the current mass, as the spin drops.
# Potential: AXION_POTENTIAL=cosine (default) or free. With the free potential the two modes are uncoupled and f_a
# enters only through the initial cloud in physical units (κ𝒩₀ with κ = (f_a/M_Pl)²); outputs get the prefix spindown_free_.
# Writes results/spindown_a<spin>_alpha0.3_fa<f_a/M_Pl>_{diagnostics,solution}.csv. The verbose log (one line per
# rhs evaluation, including b and j) is the fallback record if the run stops early.
using AxionCloud, Printf

length(ARGS) in (1, 2, 3, 4) ||
    error("usage: julia --project=. scripts/evolution_spindown.jl <f_a/M_Pl> [initial spin] [final t/M₀] [restart CSV]")
decay_constant_ratio = parse(Float64, ARGS[1])
spin_label = length(ARGS) >= 2 ? ARGS[2] : "0.9"
initial_spin = parse(Float64, spin_label)
final_time = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 1.5e9   # default 1.5e9 as in the fixed-background baseline
restart_path = length(ARGS) == 4 ? ARGS[4] : nothing

background = KerrBackground(initial_spin)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
potential_name = get(ENV, "AXION_POTENTIAL", "cosine")
potential = Dict("cosine" => CosinePotential(), "free" => FreePotential())[potential_name]
problem = RadialProblem(background, potential, mass_coupling, mode_set; mesh)
target_actions = [3e-6, 3e-6]
probe = solve_radial(problem, [1e-6, 1e-6])
amplitudes = probe.amplitudes .* sqrt.(target_actions ./ phase_actions(probe))
seeded = solve_radial(problem, amplitudes; guess=probe)
@printf("seed: A = (%.6e, %.6e), N = (%.6e, %.6e)\n", amplitudes..., phase_actions(seeded)...)
start_time, initial_background, initial_guess = 0.0, BackgroundState(1.0, initial_spin), nothing
if restart_path !== nothing
    last_row = parse.(Float64, split(readlines(restart_path)[end], ','))   # time, A1, A2, b, j
    start_time, amplitudes = last_row[1], last_row[2:3]
    initial_background = BackgroundState(last_row[4], last_row[5])
    # same problem as evolve_with_spindown's current_problem: spin j/b², α = α₀ b, mesh reused in current-mass units
    restart_problem = with_background(problem, KerrBackground(last_row[5] / last_row[4]^2), mass_coupling * last_row[4])
    for factor in exp.(range(log(1e-3), 0.0; length=30))   # amplitude continuation from the linear regime
        global initial_guess = solve_radial(restart_problem, amplitudes .* factor; guess=initial_guess)
        initial_guess.converged || error("continuation to the restart state failed at factor $factor")
    end
    @printf("restart from %s: t = %.6e, A = (%.6e, %.6e), b = %.12f, j = %.12f (a = %.6f)\n", restart_path, start_time,
            amplitudes..., last_row[4], last_row[5], last_row[5] / last_row[4]^2)
end
start_time < final_time || error("restart time $start_time is not before the final time $final_time")
span = (start_time, final_time)
prefix = joinpath(@__DIR__, "..", "results", "spindown_$(potential_name == "free" ? "free_" : "")a$(spin_label)_alpha0.3_fa$(ARGS[1])" *
                  (length(ARGS) ≥ 3 ? "_t$(ARGS[3])" : "") * (restart_path === nothing ? "" : "_restart"))
println("potential: $potential_name")
@printf("a₀ = %.4f; mesh: %d nodes, r* ∈ [%.2f, %.2f]; balance = :action; f_a/M_Pl = %.3e (f_a = %.3e GeV)\n",
        initial_spin, length(mesh.rstar), mesh.rstar[1], mesh.rstar[end], decay_constant_ratio, decay_constant_ratio * 1.2209e19)
flush(stdout)

t = @elapsed result = evolve_with_spindown(problem, amplitudes, span, decay_constant_ratio; initial_background, initial_guess,
                                           saveat=(final_time - start_time) / 60, verbose=true)
println("done: $(length(result.diagnostics)) rhs evals in $(round(t, digits=1))s"); flush(stdout)

open(prefix * "_diagnostics.csv", "w") do io
    println(io, "time,A1,A2,Re_omega1,Re_omega2,Im_omega1,Im_omega2,E,J,L_E,L_J,N1,N2,L_N1,L_N2,b,a,w")
    for d in result.diagnostics
        row = (d.time, d.amplitudes..., real.(d.mode_frequencies)..., imag.(d.mode_frequencies)..., d.charges.energy,
               d.charges.angular_momentum, d.losses.energy, d.losses.angular_momentum, d.actions..., d.action_losses...,
               d.black_hole_mass, d.black_hole_spin, d.metric_work)
        println(io, join((@sprintf("%.10e", value) for value in row), ","))
    end
end
open(prefix * "_solution.csv", "w") do io
    println(io, "time,A1,A2,b,j")
    trajectory = amplitude_trajectory(result)
    for (index, time) in enumerate(trajectory.times)
        @printf(io, "%.15e,%.15e,%.15e,%.15e,%.15e\n", time, trajectory.amplitudes[index, 1], trajectory.amplitudes[index, 2],
                trajectory.black_hole_mass[index], trajectory.black_hole_angular_momentum[index])
    end
end
println("wrote $(prefix)_{diagnostics,solution}.csv"); flush(stdout)
