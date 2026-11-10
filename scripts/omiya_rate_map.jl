# Map of where the deformation neglected by Omiya et al. 2211.01949 (Sec. III) changes the rates: joint (mutually
# deforming) action losses versus Omiya-type ones at equal occupations, on a grid of fundamental energies
# (E₂₁₁, E₃₂₂) in f_a²M, a = 0.99, cosine potential, |211⟩ + |322⟩.
#
# Omiya-type rates at occupations (N₁, N₂), as validated in scripts/omiya_flux_extrapolation.jl: solve with a tiny
# |322⟩ (E₃₂₂ ≈ 1e-2) at the A₁ whose |211⟩ phase action is N₁, then scale each Fourier component's loss by
# (N₂ / N₂,ref)^|q₂|, its leading power in the |322⟩ amplitude (their model: |211⟩ independent of A₂, |322⟩ linear on
# it, radiative modes at the order they are first sourced).
#
# Rows are appended and flushed as they finish. Continuation: columns of increasing E₂₁₁, each started from the
# previous column's reference solve; within a column, increasing E₃₂₂. A column stops at its first failed point.
#
# Run from the package root:  julia --project=. scripts/omiya_rate_map.jl <α> [quick | column]   (submit to Slurm;
# `quick` runs one point as a smoke test; an integer `column` computes only that E₂₁₁ column, after stepping the
# reference through the lower columns, and writes `..._col<column>.csv`, so the columns can run as parallel jobs)
using AxionCloud, Printf

mass_coupling = parse(Float64, ARGS[1])
quick = length(ARGS) ≥ 2 && ARGS[2] == "quick"
only_column = length(ARGS) ≥ 2 && !quick ? parse(Int, ARGS[2]) : nothing
primary_energies = quick ? [300.0] : [30.0, 100.0, 300.0, 1000.0, 2000.0, 2700.0]
secondary_energies = quick ? [10.0] : [1.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]
reference_secondary_energy = 1e-2

background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, CosinePotential(), mass_coupling, mode_set; mesh)
components = mode_set.components
secondary_power = [abs(component.q[2]) for component in components]
weights = AxionCloud.phase_weights(mode_set)
output_path = joinpath(@__DIR__, "..", "results",
                       @sprintf("omiya_rate_map_a0.99_alpha%.2f%s.csv", mass_coupling,
                                quick ? "_quick" : only_column === nothing ? "" : "_col$(only_column)"))
@printf("α = %.2f, mesh: %d nodes, r* ∈ [%.2f, %.2f]; output %s\n", mass_coupling, length(mesh.rstar),
        mesh.rstar[1], mesh.rstar[end], output_path)
flush(stdout)

"""Outward action flux of each Fourier component (unweighted by phase) at mesh `node`; same lines as `action_flux`."""
function component_action_fluxes(solution, node)
    state = AxionCloud.physical_state(solution.layout, solution.unknowns, solution.scale, node)
    wronskian = imag.(conj.(state.profiles) .* state.slopes)
    prefactor = 4pi * AxionCloud.r2a2(problem.background, problem.mesh.offsets[node])
    [prefactor * sum(wronskian[channels]) for channels in mode_set.component_channels]
end
component_losses(solution) = component_action_fluxes(solution, solution.layout.n_nodes) -
                             component_action_fluxes(solution, 1)

"""
Small-secondary solve whose |211⟩ phase action is `primary_action`: secant iteration in ln A₁ with the measured slope
d ln N₁ / d ln A₁ (2 for a linear cloud; it falls towards 0 near the end of the |211⟩ branch, where the fixed point
A₁ ← A₁ √(N₁/N₁(A₁)) stalls: α = 0.42 at E₂₁₁ = 1000, slope ≈ 0.8), each step limited to a factor 1.5 in A₁.
"""
function primary_matched_reference(primary_action, start)
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

labels = ["q$(component.q[1])_$(component.q[2])" for component in components]
open(output_path, "w") do io
    println(io, join(["target_E211", "target_E322", "E211", "E322", "Re_omega1", "Re_omega2", "N1", "N2",
                      "joint_L_N1", "joint_L_N2", "omiya_L_N1", "omiya_L_N2",
                      ("joint_loss_" .* labels)..., ("omiya_loss_" .* labels)...], ","))
end

column_start = nothing
for (column_index, primary_energy) in enumerate(primary_energies)
    only_column !== nothing && column_index > only_column && break
    # Column reference: |211⟩ at the target energy with a tiny |322⟩; it starts both the joint continuation and
    # the Omiya-type matching.
    reference = try
        column_start === nothing ? amplitudes_for_energy(problem, [primary_energy, reference_secondary_energy]) :
            amplitudes_for_energy(problem, [primary_energy, reference_secondary_energy];
                                  initial=column_start.amplitudes, guess=column_start)
    catch error_value
        @printf("E₂₁₁ = %g: reference calibration failed (%s); stopping\n", primary_energy, sprint(showerror, error_value))
        break
    end
    if !reference.converged
        @printf("E₂₁₁ = %g: reference calibration did not converge; stopping\n", primary_energy)
        break
    end
    global column_start = reference
    @printf("column E₂₁₁ = %g: reference A = (%.6e, %.6e)\n", primary_energy, reference.amplitudes...)
    flush(stdout)
    only_column !== nothing && column_index < only_column && continue

    previous = reference
    for secondary_energy in secondary_energies
        elapsed = @elapsed row = try
            joint = amplitudes_for_energy(problem, [primary_energy, secondary_energy]; initial=[previous.amplitudes[1],
                                          previous.amplitudes[2] * sqrt(secondary_energy / fundamental_energies(previous)[2])],
                                          guess=previous)
            joint.converged || error("joint calibration did not converge")
            actions = phase_actions(joint)
            matched = primary_matched_reference(actions[1], reference)
            action_ratio = actions[2] / phase_actions(matched)[2]
            omiya_components = component_losses(matched) .* action_ratio .^ secondary_power
            (; joint, actions, joint_components=component_losses(joint), omiya_components,
             omiya_losses=weights * omiya_components)
        catch error_value
            @printf("  E₃₂₂ = %g: failed (%s); next column\n", secondary_energy, sprint(showerror, error_value))
            nothing
        end
        row === nothing && break
        previous = row.joint
        joint_losses = action_losses(row.joint)
        values = (primary_energy, secondary_energy, fundamental_energies(row.joint)..., real.(mode_frequencies(row.joint))...,
                  row.actions..., joint_losses..., row.omiya_losses..., row.joint_components..., row.omiya_components...)
        open(output_path, "a") do io
            println(io, join((@sprintf("%.10e", value) for value in values), ","))
        end
        @printf("  E₃₂₂ = %-6g N = (%.4e, %.4e)  joint/Omiya − 1 = (%+.3e, %+.3e)  [%.0f s]\n", secondary_energy,
                row.actions..., (joint_losses ./ row.omiya_losses .- 1)..., elapsed)
        flush(stdout)
    end
end
println("done")
