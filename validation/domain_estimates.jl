# Print the physically motivated domains at a few parameter points, and compare with the
# fixed defaults. Run from the package root:  julia --project=. validation/domain_estimates.jl
using AxionCloud, Printf
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
for (spin, mass_coupling) in ((0.9, 0.3), (0.99, 0.3), (0.99, 0.2), (0.9, 0.4), (0.999, 0.3))
    background = KerrBackground(spin)
    estimate = domain_estimate(background, mass_coupling, mode_set)
    mesh = RadialMesh(background, mass_coupling, mode_set)
    @printf("a = %5.3f α = %.2f: r* ∈ [%7.1f, %8.1f]  containment r = %7.1f  normalization limit r = %7.1f  asymptotic r = %7.1f%s  κ_min = %.4f  k_min = %.4f  peak r = %5.1f  nodes = %d\n",
            spin, mass_coupling, estimate.rstar_min, estimate.rstar_max, estimate.containment_radius, estimate.normalization_radius,
            estimate.asymptotic_radius, estimate.capped ? " (capped)" : "", estimate.slowest_decay,
            estimate.smallest_wavenumber, estimate.peak_radius, length(mesh.rstar))
end
