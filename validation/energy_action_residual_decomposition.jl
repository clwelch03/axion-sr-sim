# Where does D_k = ∂_{A_k}E − Σ_i ω_i ∂_{A_k}N_i come from? (reviewer item H, docs/loose-ends.txt)
#
# The code's E is the Killing energy −∫T^t_t √−g, which for a scalar field equals the canonical energy
# Σ_i ω_i N_i − L̄ term by term, with L̄ the phase-averaged Lagrangian on the same mesh and quadrature.
# Differentiating along the solution family at fixed background and integrating the radial kinetic term
# by parts gives, with physical profiles R, tortoise slopes Q = dR/dr*, and real component frequencies ω_q,
#
#     D_k = S_k + G_k,
#     S_k = 4π [ (r²+a²) Σ_{l,q} Re(Q* ∂_{A_k}R) ]_{r_in}^{r_out}                       (surface term)
#     G_k = 4π Re ∫ dr* (r²+a²) Σ_{l,q} δQ' ∂_{A_k}R*                                  (growth term)
#
# where δQ' = (radial right-hand side dQ/dr* at the real ω_q) − (the same at the complex ω̃_q the profiles
# solve). S_k vanishes only if the endpoints carry no field; G_k vanishes only for real frequencies. The
# remainder D_k − S_k − G_k is quadrature error of the discretised integrals (it should shrink with the step).
#
# Usage:  julia --project=. validation/energy_action_residual_decomposition.jl <point 1|2|3> <mesh 1|2|3>
#   points: 1 = A (1e-6, 1e-7) linear; 2 = (1e-2, 1e-4); 3 = (1e-2, 1e-3)
#   meshes: 1 = r* ∈ [−60, 500], step 0.4 (as validation/domain_actions_check.jl); 2 = same, step 0.2;
#           3 = r* ∈ [−60, 700], step 0.4
# a = 0.9, α = 0.3, modes |211⟩ and |322⟩, default (conjugate) driven-growth rule. Potential from AXION_POTENTIAL
# (cosine, the default, or chiral).
using AxionCloud, Printf, LinearAlgebra
const AC = AxionCloud

point_index = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
mesh_index = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
amplitude_targets = ([1e-6, 1e-7], [1e-2, 1e-4], [1e-2, 1e-3])
mesh_settings = ((rstar_max=500.0, step=0.4), (rstar_max=500.0, step=0.2), (rstar_max=700.0, step=0.4))

background = KerrBackground(0.9)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
settings = mesh_settings[mesh_index]
mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=settings.rstar_max, step=settings.step,
                  stretch_from=150.0, stretch_factor=1.03, max_step=8settings.step)
potential_name = get(ENV, "AXION_POTENTIAL", "cosine")
potential = Dict("cosine" => CosinePotential(), "chiral" => ChiralPotential())[potential_name]
problem = RadialProblem(background, potential, mass_coupling, mode_set; mesh)
target = amplitude_targets[point_index]

# amplitude continuation from the linear regime, as in validation/domain_actions_check.jl
solve_time = @elapsed begin
    global solution = solve_radial(problem, [1e-6, 1e-7])
    if maximum(target) >= 1e-4
        for fraction in exp10.(range(-4, 0; length=13)[2:end])
            global solution = solve_radial(problem, target .* fraction; guess=solution)
            solution.converged || @warn "non-converged continuation step" fraction solution.residual_norm
        end
    end
end
layout, unknowns, scale = solution.layout, solution.unknowns, solution.scale
@printf("%s; point %d, A = (%.1e, %.1e); mesh %d: r* ∈ [%.0f, %.0f], step %.2f, %d nodes; converged %s (%.0f s)\n",
        potential_name, point_index, target..., mesh_index, mesh.rstar[1], mesh.rstar[end], settings.step, layout.n_nodes,
        solution.converged, solve_time)

# D_k from the code's own charge and action Jacobians
directions = amplitude_response(solution)
energy_wrt_amplitudes = charge_jacobian(solution; directions)[1, :]
actions_wrt_amplitudes = action_jacobian(solution; directions)           # [mode i, amplitude k]
real_mode_frequencies = real.(mode_frequencies(solution))
residual = energy_wrt_amplitudes .- vec(real_mode_frequencies' * actions_wrt_amplitudes)

# surface and growth terms
complex_component_frequencies = AC.component_frequencies(mode_set, mode_frequencies(solution), problem.driven_growth_rule)
real_component_frequencies = real.(complex_component_frequencies)
profile_derivative(node, amplitude_index) = scale .* AC.NodeState(directions[AC.node_range(layout, node), amplitude_index]).profiles

function surface_density(node, amplitude_index)
    state = AC.physical_state(layout, unknowns, scale, node)
    return 4pi * AC.r2a2(background, mesh.offsets[node]) *
           sum(real.(conj.(state.slopes) .* profile_derivative(node, amplitude_index)))
end

function growth_density(node, amplitude_index)
    offset = mesh.offsets[node]
    state = AC.NodeState(unknowns[AC.node_range(layout, node)])   # scaled profiles, as radial_rhs expects
    _, slope_rate_real = AC.radial_rhs(problem, offset, state.profiles, state.slopes, real_component_frequencies, scale)
    _, slope_rate_complex = AC.radial_rhs(problem, offset, state.profiles, state.slopes, complex_component_frequencies, scale)
    slope_rate_difference = scale .* (slope_rate_real .- slope_rate_complex)
    return 4pi * AC.r2a2(background, offset) *
           sum(real.(slope_rate_difference .* conj.(profile_derivative(node, amplitude_index))))
end

simpson = AC.simpson_weights(mesh.rstar)
for amplitude_index in 1:2
    surface_outer = surface_density(layout.n_nodes, amplitude_index)
    surface_inner = surface_density(1, amplitude_index)
    surface = surface_outer - surface_inner
    growth = sum(simpson[node] * growth_density(node, amplitude_index) for node in 1:layout.n_nodes)
    scale_E = abs(energy_wrt_amplitudes[amplitude_index])
    @printf("  k = %d:  ∂E/∂A = %+.6e   D = %+.6e (D/|∂E| = %+.2e)\n", amplitude_index,
            energy_wrt_amplitudes[amplitude_index], residual[amplitude_index], residual[amplitude_index] / scale_E)
    @printf("         S_out = %+.6e   S_in = %+.6e   S = %+.6e   G = %+.6e\n", surface_outer, surface_inner, surface, growth)
    @printf("         D − S − G = %+.3e   (relative to |D|: %+.2e; to |∂E|: %+.2e)\n",
            residual[amplitude_index] - surface - growth,
            (residual[amplitude_index] - surface - growth) / abs(residual[amplitude_index]),
            (residual[amplitude_index] - surface - growth) / scale_E)
end
@printf("  ω = (%.10f, %.10f), γ = (%+.4e, %+.4e)\n", real_mode_frequencies..., imag.(mode_frequencies(solution))...)
