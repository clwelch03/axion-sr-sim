# Quasi-equilibrium (plateau) of the joint |211⟩ + |322⟩ cloud on a fixed Kerr background, found directly instead of
# by evolution, and the error of Omiya-type rates there (Omiya et al. 2211.01949 Sec. III: |211⟩ independent of |322⟩,
# |322⟩ linear on it; built as in scripts/omiya_rate_map.jl), over a sweep of spins at one α.
#
# Plateau: amplitudes A = (A₁, A₂) with both action losses zero. Unknowns x = ln A_i; residual r_i = L_N,i / N_i (the
# rate of ln N_i), which excludes the trivial root A₂ = 0. Newton with the exact Jacobian ∂L/∂A of the coupled solution
# (`observable_jacobian`) and ∂N/∂A (`action_jacobian`); steps limited to 0.5 in each ln A_i.
#
# Error at the plateau N*:
#   measure      = max_i |L_N,i^Omiya(N*)| / G_i, with G_i = Σ_q |q_i| |loss_q| the gross action flux through phase i
#                  (the rate map's colour, here with L_joint = 0);
#   plateau shift δN = −(∂L/∂N)⁻¹ L^Omiya(N*), first order, with the joint ∂L/∂N = (∂L/∂A)(∂N/∂A)⁻¹: an estimate of
#                  where Omiya's plateau sits relative to ours.
#
# Run from the package root:  julia --project=. scripts/plateau_scan.jl <α> [spins...]   (default spins 0.99 0.95 0.9
# 0.85 0.8; each spin starts from the previous plateau's amplitudes)
using AxionCloud, LinearAlgebra, Printf

mass_coupling = parse(Float64, ARGS[1])
spins = length(ARGS) ≥ 2 ? parse.(Float64, ARGS[2:end]) : [0.99, 0.95, 0.9, 0.85, 0.8]
alpha_label = round(mass_coupling; digits=2) == mass_coupling ? @sprintf("%.2f", mass_coupling) : @sprintf("%.3f", mass_coupling)
output_path = joinpath(@__DIR__, "..", "results", "plateau_scan_alpha$(alpha_label).csv")
# inner_energy_flux, inner_angular_momentum_flux: outward fluxes at the innermost node (the horizon side), in f_a²
# units of the current mass (`outward_flux(solution, 1)`); positive means the hole loses E or J (added 2026-11-10,
# absent from older CSVs).
open(output_path, "w") do io
    println(io, "alpha,spin,converged,iterations,A1,A2,N1,N2,E211,E322,Re_omega1,Re_omega2,growth_rate_211," *
                "gross_1,gross_2,omiya_L_N1,omiya_L_N2,measure,shift_N1_rel,shift_N2_rel," *
                "inner_energy_flux,inner_angular_momentum_flux")
end

"""Outward action flux of each Fourier component at mesh `node`, generic in the unknowns (for ForwardDiff)."""
function component_action_fluxes(problem, layout, unknowns, scale, node)
    state = AxionCloud.physical_state(layout, unknowns, scale, node)
    wronskian = imag.(conj.(state.profiles) .* state.slopes)
    prefactor = 4pi * AxionCloud.r2a2(problem.background, problem.mesh.offsets[node])
    [prefactor * sum(wronskian[channels]) for channels in problem.mode_set.component_channels]
end
component_losses(problem, layout, unknowns, scale) =
    component_action_fluxes(problem, layout, unknowns, scale, layout.n_nodes) -
    component_action_fluxes(problem, layout, unknowns, scale, 1)
component_losses(solution) = component_losses(solution.problem, solution.layout, solution.unknowns, solution.scale)
phase_losses(problem, layout, unknowns, scale) =
    AxionCloud.phase_weights(problem.mode_set) * component_losses(problem, layout, unknowns, scale)

"""Newton solve for the plateau from the converged solution `start`; returns (solution, converged, iterations)."""
function find_plateau(problem, start::RadialSolution; maxiter=40, tol=1e-7)
    solution = start
    for iteration in 1:maxiter
        solution.converged || return solution, false, iteration
        actions, losses = phase_actions(solution), action_losses(solution)
        weights = AxionCloud.phase_weights(problem.mode_set)
        gross = [sum(abs.(weights[phase, :]) .* abs.(component_losses(solution))) for phase in 1:2]
        maximum(abs.(losses) ./ gross) < tol && return solution, true, iteration
        directions = amplitude_response(solution)
        loss_jacobian = AxionCloud.observable_jacobian(solution, phase_losses; directions)          # ∂L_i/∂A_j
        action_derivatives = action_jacobian(solution; directions)                      # ∂N_i/∂A_j
        residual = losses ./ actions
        jacobian = (loss_jacobian ./ actions .- losses ./ actions .^ 2 .* action_derivatives) .* solution.amplitudes'
        step = clamp.(-(jacobian \ residual), -0.5, 0.5)
        new_amplitudes = solution.amplitudes .* exp.(step)
        new_solution = solve_radial(problem, new_amplitudes; guess=solution)
        backoff = 0
        while !new_solution.converged && backoff < 4
            backoff += 1
            step ./= 2
            new_amplitudes = solution.amplitudes .* exp.(step)
            new_solution = solve_radial(problem, new_amplitudes; guess=solution)
        end
        solution = new_solution
        @printf("    Newton %2d: A = (%.6e, %.6e)  L/G = (%+.2e, %+.2e)\n", iteration, solution.amplitudes...,
                (losses ./ gross)...)
        flush(stdout)
    end
    return solution, false, maxiter
end

"""Small-secondary solve whose |211⟩ phase action is `primary_action` (secant in ln A₁, as in omiya_rate_map.jl)."""
function primary_matched_reference(problem, primary_action, start)
    solution = start
    primary, secondary = start.amplitudes
    log_slope = 2.0
    for _ in 1:40
        current_action = phase_actions(solution)[1]
        abs(current_action / primary_action - 1) < 1e-9 && return solution
        step = clamp(log(primary_action / current_action) / log_slope, -log(1.5), log(1.5))
        new_primary = primary * exp(step)
        new_solution = solve_radial(problem, [new_primary, secondary]; guess=solution)
        new_solution.converged || error("matched reference solve did not converge")
        abs(step) > 1e-12 && (log_slope = clamp(log(phase_actions(new_solution)[1] / current_action) / step, 0.02, 3.0))
        primary, solution = new_primary, new_solution
    end
    error("A₁ matching did not converge for N₁ = $primary_action")
end

mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
start_amplitudes = nothing
start_energies = [50.0, 80.0]
for spin in spins
    background = KerrBackground(spin)
    mesh = RadialMesh(background, mass_coupling, mode_set)
    problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
    growth_rate = growth_rates(solve_radial(with_potential(problem, FreePotential()), [1e-6, 1e-6]))[1]
    @printf("α = %.2f, a = %.3f: mesh %d nodes, r* ∈ [%.2f, %.2f]; linear growth rates ω₂₁₁,I = %+.4e\n",
            mass_coupling, spin, length(mesh.rstar), mesh.rstar[1], mesh.rstar[end], growth_rate)
    flush(stdout)
    if growth_rate ≤ 0
        @printf("  a = %.3f: |211⟩ not superradiant, no two-cloud plateau sought\n", spin)
        open(output_path, "a") do io
            println(io, join((@sprintf("%.4f", mass_coupling), @sprintf("%.4f", spin), "-1", fill("NaN", 19)...), ","))
        end
        continue
    end
    elapsed = @elapsed result = try
        # Start: a cold solve at the previous spin's plateau amplitudes, kept only if it converged to energies within
        # a factor 2 of the previous plateau's (at α = 0.187 a cold solve at E = (5, 8) converged towards a spurious
        # state with E₂₁₁ ≈ 1e-7). Otherwise (or at the first spin) a two-stage calibration: E = (0.5, 0.8) from the
        # linear estimate, then the previous plateau energies (default 50, 80) warm-started from it. At α = 0.187 the
        # one-stage calibration to (50, 80) from the linear estimate does not converge; the two-stage one does.
        start = start_amplitudes === nothing ? nothing : solve_radial(problem, start_amplitudes)
        if start !== nothing && start.converged
            energy_ratios = fundamental_energies(start) ./ start_energies
            all(0.5 .< energy_ratios .< 2) || (start = nothing)
        end
        if start === nothing || !start.converged
            small = amplitudes_for_energy(problem, [0.5, 0.8])
            small.converged || error("small-energy start calibration did not converge")
            start = amplitudes_for_energy(problem, start_energies; initial=small.amplitudes, guess=small)
            start.converged || error("start calibration did not converge")
        end
        solution, converged, iterations = find_plateau(problem, start)
        converged || error("plateau Newton did not converge in $iterations iterations")
        actions = phase_actions(solution)
        weights = AxionCloud.phase_weights(mode_set)
        gross = [sum(abs.(weights[phase, :]) .* abs.(component_losses(solution))) for phase in 1:2]
        # Omiya-type rates at the plateau occupations.
        energies = real.(mode_frequencies(solution)) .* actions
        small_secondary = solution.amplitudes[2] * sqrt(1e-2 / energies[2])
        reference = solve_radial(problem, [solution.amplitudes[1], small_secondary]; guess=solution)
        matched = primary_matched_reference(problem, actions[1], reference)
        action_ratio = actions[2] / phase_actions(matched)[2]
        omiya_losses = weights * (component_losses(matched) .* action_ratio .^ [abs(c.q[2]) for c in mode_set.components])
        directions = amplitude_response(solution)
        loss_jacobian = AxionCloud.observable_jacobian(solution, phase_losses; directions)
        action_derivatives = action_jacobian(solution; directions)
        shift = -((loss_jacobian / action_derivatives) \ omiya_losses)
        (; solution, iterations, actions, energies, gross, omiya_losses, measure=maximum(abs.(omiya_losses) ./ gross),
         shift_relative=shift ./ actions)
    catch error_value
        @printf("  a = %.3f: failed (%s)\n", spin, first(sprint(showerror, error_value), 300))
        nothing
    end
    if result === nothing
        open(output_path, "a") do io
            println(io, join((@sprintf("%.4f", mass_coupling), @sprintf("%.4f", spin), "0", fill("NaN", 19)...), ","))
        end
        continue
    end
    global start_amplitudes = result.solution.amplitudes
    global start_energies = fundamental_energies(result.solution)
    inner_flux = outward_flux(result.solution, 1)
    values = (mass_coupling, spin, 1, result.iterations, result.solution.amplitudes..., result.actions...,
              result.energies..., real.(mode_frequencies(result.solution))..., growth_rate, result.gross...,
              result.omiya_losses..., result.measure, result.shift_relative...,
              inner_flux.energy, inner_flux.angular_momentum)
    open(output_path, "a") do io
        println(io, join((@sprintf("%.10e", value) for value in values), ","))
    end
    @printf("  a = %.3f: plateau ωN = (%.4g, %.4g) after %d Newton steps; measure = %.3e; Omiya plateau shift δN/N = (%+.3e, %+.3e)  [%.0f s]\n",
            spin, result.energies..., result.iterations, result.measure, result.shift_relative..., elapsed)
    flush(stdout)
end
println("done")
