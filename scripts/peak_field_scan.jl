# Peak field |Φ|/f_a of a single |211⟩ cloud along its nonlinear branch, against the cloud energy, at a = 0.99,
# α = 0.3. The field is Φ = Σ_q F_q e^{-iq·ξ} + c.c. in units of f_a (see `reconstruct_field`); the peak is the maximum
# over the radial mesh, a 13-node Gauss–Legendre polar grid (odd, so it contains the equator) and the phase grid.
#
# Continuation in the asymptotic amplitude A (geometric steps, each solve warm-started from the previous one), not in
# E, so the scan passes the branch end where dE/dA → 0. The step is halved when Newton fails, and the scan stops
# after the step falls below 1e-3 in ln A.
#
# Run from the package root:
#   julia --project=scripts scripts/peak_field_scan.jl <cosine|chiral> [A_start] [A_stop]
# Writes results/peak_field_scan_<potential>.csv.
using AxionCloud, Printf

potential_name = ARGS[1]
amplitude_start = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 3e-3
amplitude_stop = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.5
potential = Dict("cosine" => CosinePotential(), "chiral" => ChiralPotential())[potential_name]
spin, mass_coupling = 0.99, 0.3

background = KerrBackground(spin)
mode_set = ModeSet([BoundMode(2, 1, 1)])
mesh = RadialMesh(background, mass_coupling, mode_set)
problem = RadialProblem(background, potential, mass_coupling, mode_set; mesh)
equator_quadrature = FieldQuadrature(mode_set; n_polar=13)

"""Maximum of |Φ|/f_a over the mesh, the polar grid and the phases, and the tortoise coordinate where it occurs."""
function peak_field(solution)
    node_profiles = profiles(solution)
    peaks = [maximum(abs, reconstruct_field(equator_quadrature, node_profiles[node, :])) for node in axes(node_profiles, 1)]
    node = argmax(peaks)
    return peaks[node], mesh.rstar[node]
end

output_path = joinpath(@__DIR__, "..", "results", "peak_field_scan_$(potential_name).csv")
open(output_path, "w") do io
    println(io, "# single |211⟩, a = $spin, α = $mass_coupling, $potential_name potential, mesh $(length(mesh.rstar)) nodes")
    println(io, "amplitude,fundamental_energy,total_energy,angular_momentum,Re_omega,growth_rate,peak_field,peak_rstar,residual")
end
@printf("%s: mesh %d nodes, r* ∈ [%.1f, %.1f]\n", potential_name, length(mesh.rstar), mesh.rstar[1], mesh.rstar[end])

amplitude, log_step, previous = amplitude_start, log(1.1), nothing
while amplitude ≤ amplitude_stop
    started = time()
    solution = solve_radial(problem, [amplitude]; guess=previous)
    if !solution.converged
        global log_step /= 2
        log_step < 1e-3 && (println("stopping: no convergence at A = $amplitude"); break)
        global amplitude = previous.amplitudes[1] * exp(log_step)
        continue
    end
    energy = fundamental_energies(solution)[1]
    charges = cloud_charges(solution)
    frequency = mode_frequencies(solution)[1]
    field, field_rstar = peak_field(solution)
    @printf("A = %.5f  E_fund = %9.3f  E = %9.3f  ω = %.8f  Γ = %.3e  |Φ|max/f_a = %.4f at r* = %.1f  (%.0f s)\n",
            amplitude, energy, charges.energy, real(frequency), growth_rates(solution)[1], field, field_rstar,
            time() - started)
    flush(stdout)
    open(output_path, "a") do io
        @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.3g\n", amplitude, energy, charges.energy,
                charges.angular_momentum, real(frequency), growth_rates(solution)[1], field, field_rstar,
                solution.residual_norm)
    end
    global previous = solution
    global amplitude *= exp(log_step)
end
println(output_path)
