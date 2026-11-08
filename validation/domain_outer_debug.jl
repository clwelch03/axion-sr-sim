# Verbose look at the joint solve on a far-out stretched mesh (why was r*_max = 1000 slow?).
# Run from the package root:  julia --project=. validation/domain_outer_debug.jl
using AxionCloud, Printf
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
rstar_max = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 1000.0
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step=3.7)
println("$(length(mesh.rstar)) nodes, largest step $(maximum(diff(mesh.rstar)))"); flush(stdout)
free = RadialProblem(background, FreePotential(), 0.3, mode_set; mesh)
t = @elapsed linear = solve_radial(free, [1e-6, 1e-6]; verbose=true)
println("linear: $(round(t, digits=1))s converged=$(linear.converged) ω = ", mode_frequencies(linear)); flush(stdout)
estimate = linear.amplitudes .* sqrt.([650.0, 650.0] ./ fundamental_energies(linear))
println("linear estimate A = ", estimate); flush(stdout)
cosine = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
t = @elapsed joint = solve_radial(cosine, estimate; guess=linear, verbose=true)
println("first nonlinear solve: $(round(t, digits=1))s converged=$(joint.converged) iterations=$(joint.iterations) ω = ", mode_frequencies(joint)); flush(stdout)
t = @elapsed calibrated = amplitudes_for_energy(cosine, [650.0, 650.0]; initial=joint.amplitudes, guess=joint, verbose=true)
println("calibration: $(round(t, digits=1))s ω = ", mode_frequencies(calibrated), " L = ", boundary_losses(calibrated)); flush(stdout)
