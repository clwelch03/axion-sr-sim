# Per-component diagnostics of the joint a = 0.99 solution: frequencies, wavenumbers, outer
# coefficients and fluxes, which identify the near-threshold (−1, 2) emission and the m = 0
# horizon absorption. Run from the package root:  julia --project=. validation/component_diagnostics.jl
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=200.0, step=0.4)
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
joint = amplitudes_for_energy(problem, [650.0, 650.0])
println("ω = ", mode_frequencies(joint))
frequencies = AC.component_frequencies(mode_set, mode_frequencies(joint))
R, Q = profiles(joint), slopes(joint)
outer_prefactor = 4pi * AC.r2a2(background, mesh.offsets[end])
inner_prefactor = 4pi * AC.r2a2(background, mesh.offsets[1])
for (index, component) in enumerate(mode_set.components)
    frequency = frequencies[index]
    k = outer_wavenumber(frequency, 0.3)
    channel = mode_set.component_channels[index][1]
    outer_flux = outer_prefactor * real(frequency) * imag(conj(R[end, channel]) * Q[end, channel])
    inner_flux = inner_prefactor * real(frequency) * imag(conj(R[1, channel]) * Q[1, channel])
    @printf("q = %-8s m = %2d  ω_q = %.6f%+.2ei  k = %.4f%+.4fi  ν = %s  |R_end| = %.2e  F_E(out) = %+.3e  F_E(in) = %+.3e  peak|R| = %.3e\n",
            string(component.q), component.m, real(frequency), imag(frequency), real(k), imag(k),
            round(AC.outer_exponent(frequency, 0.3, k); digits=2), abs(R[end, channel]), outer_flux, inner_flux, maximum(abs, R[:, channel]))
end
println("boundary losses: ", boundary_losses(joint))
