# Sensitivity of the joint a = 0.99, α = 0.3, E₁ = E₂ = 650 rates to the series orders (J_H, J_B),
# the radial step and the quadrature sizes, on the stretched r* ∈ [−100, 400] mesh.
# Run from the package root:  julia --project=. validation/truncation_convergence.jl  (≈ 1 h)
using AxionCloud, Printf
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
mesh_for(step) = RadialMesh(background; rstar_min=-100.0, rstar_max=400.0, step, stretch_from=150.0, stretch_factor=1.03, max_step=3.7)
reference = nothing
function run(label; step=0.4, horizon_order=3, outer_order=7, n_polar=12, n_phase=16)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh=mesh_for(step), horizon_order, outer_order, n_polar, n_phase)
    t = @elapsed joint = amplitudes_for_energy(problem, [650.0, 650.0])
    rates = imag.(mode_frequencies(joint))
    global reference
    reference === nothing && (reference = rates)
    @printf("%-34s %4.0fs  Im ω = (%.6e, %.6e)  rel. change vs reference = (%+.2e, %+.2e)  A = %s\n",
            label, t, rates[1], rates[2], rates[1] / reference[1] - 1, rates[2] / reference[2] - 1, joint.amplitudes)
    flush(stdout)
end
run("reference (J_H=3, J_B=7, h=0.4, 12/16)")
run("J_H = 5, J_B = 9"; horizon_order=5, outer_order=9)
run("J_H = 7, J_B = 11"; horizon_order=7, outer_order=11)
run("step 0.2"; step=0.2)
run("n_polar = 16, n_phase = 24"; n_polar=16, n_phase=24)
