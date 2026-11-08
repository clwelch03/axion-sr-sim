# Convergence of the joint a = 0.99, α = 0.3, E₁ = E₂ = 650 solution with the outer endpoint on
# stretched meshes, and of the linear two-mode rates with the inner endpoint.
# Run from the package root:  julia --project=. validation/domain_convergence.jl  (≈ 1 h)
using AxionCloud, Printf
const AC = AxionCloud
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])

println("== inner endpoint, linear two-mode problem, a = 0.9, α = 0.3 (r*_max = 180, step 0.4)")
background = KerrBackground(0.9)
for rstar_min in (-30.0, -45.0, -60.0, -80.0, -120.0)
    mesh = RadialMesh(background; rstar_min, rstar_max=180.0, step=0.4)
    problem = RadialProblem(background, FreePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    @printf("r*_min = %6.1f  s_min = %.2e  s_min/gap = %.2e  ω = %s  iterations = %d\n", rstar_min, mesh.offsets[1],
            mesh.offsets[1] / background.horizon_gap, mode_frequencies(solution), solution.iterations)
    flush(stdout)
end

println("== outer endpoint, joint cosine solution, a = 0.99, α = 0.3, E₁ = E₂ = 650")
background = KerrBackground(0.99)
for (rstar_max, max_step) in ((400.0, 3.7), (1000.0, 3.7), (2000.0, 3.7), (4000.0, 3.7))
    mesh = RadialMesh(background; rstar_min=-100.0, rstar_max, step=0.4, stretch_from=150.0, stretch_factor=1.03, max_step)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    t = @elapsed joint = amplitudes_for_energy(problem, [650.0, 650.0])
    frequencies = AC.component_frequencies(mode_set, mode_frequencies(joint))
    R, Q = profiles(joint), slopes(joint)
    index = findfirst(component -> component.q == [-1, 2], mode_set.components)
    channel = mode_set.component_channels[index][1]
    emitted = 4pi * AC.r2a2(background, mesh.offsets[end]) * real(frequencies[index]) * imag(conj(R[end, channel]) * Q[end, channel])
    @printf("r*_max = %5.0f (%4d nodes, %4.0fs, %2d it): A = %s\n   ω = %s\n   L = %s  F_E(−1,2) = %.4e  k(−1,2) = %.5f  |R322_end|/peak = %.2e\n",
            rstar_max, length(mesh.rstar), t, joint.iterations, joint.amplitudes, mode_frequencies(joint),
            collect(boundary_losses(joint)), emitted, real(outer_wavenumber(frequencies[index], 0.3)),
            abs(R[end, 2]) / maximum(abs, R[:, 2]))
    flush(stdout)
end
