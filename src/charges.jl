"""
Phase-averaged conserved charges of the interacting field, the fluxes that
carry them through the numerical boundaries, and the response of the charges
to the prescribed amplitudes.

Energies are in units of `f_a² M` and angular momenta in `f_a² M²`; fluxes are
per unit `t/M`. With the Killing vectors `∂_t`, `∂_φ` of Kerr,

    E = −∫ T^t_t √−g d³x,    J = ∫ T^t_φ √−g d³x,

and the outward fluxes through a sphere are `F_E = −∮ T^r_t √−g`, `F_J = ∮ T^r_φ √−g`,
so that `dE/dt = −(F_E(out) − F_E(in))` and likewise for `J`.
"""

"""
    Charges(energy, angular_momentum)

An `(E, J)` pair: cloud charges, outward fluxes, or boundary losses depending
on which function produced it.
"""
struct Charges{T<:Number}
    energy::T
    angular_momentum::T
end

Base.:-(left::Charges, right::Charges) = Charges(left.energy - right.energy, left.angular_momentum - right.angular_momentum)
Base.:+(left::Charges, right::Charges) = Charges(left.energy + right.energy, left.angular_momentum + right.angular_momentum)
Base.:*(factor::Number, charges::Charges) = Charges(factor * charges.energy, factor * charges.angular_momentum)
Base.collect(charges::Charges) = [charges.energy, charges.angular_momentum]

"""
    FieldBilinears

Phase-averaged quadratic quantities of the field at one radial node, each a
vector over the polar quadrature nodes. With `Φ = Σ_q F_q e^{-i q·ξ} + c.c.`
every bilinear is `Σ_q 2 Re[(…F_q)(…F_q)*]` at leading adiabatic order.
"""
struct FieldBilinears{T<:Real}
    dt_dt::Vector{T}          # ⟨(∂_t Φ)²⟩
    dphi_dphi::Vector{T}      # ⟨(∂_φ Φ)²⟩
    dr_dr::Vector{T}          # ⟨(∂_r Φ)²⟩
    dtheta_dtheta::Vector{T}  # ⟨(∂_θ Φ)²⟩
    dt_dphi::Vector{T}        # ⟨∂_t Φ ∂_φ Φ⟩
    potential::Vector{T}      # ⟨U(Φ)⟩
end

"""
    field_bilinears(problem, s, profiles, slopes, real_frequencies)

Bilinears at offset `s` from the *physical* profiles and tortoise slopes of all
channels and the real component frequencies `ω_q`.
"""
function field_bilinears(problem::RadialProblem, s, profiles::AbstractVector, slopes::AbstractVector,
                         real_frequencies::AbstractVector; potential::AxionPotential=problem.potential)
    quadrature, mode_set = problem.quadrature, problem.mode_set
    harmonics, harmonics_dtheta = quadrature.angular.harmonics, quadrature.angular.harmonics_dtheta
    azimuthal = [component.m for component in mode_set.components]

    component_values = component_profiles(quadrature, profiles)
    component_dr = component_profiles(quadrature, slopes) ./ dr_drstar(problem.background, s)
    component_dtheta = zeros(eltype(component_values), size(component_values))
    for (component_index, channels) in enumerate(mode_set.component_channels)
        component_dtheta[:, component_index] = harmonics_dtheta[:, channels] * profiles[channels]
    end

    squared = 2 .* abs2.(component_values)
    FieldBilinears(squared * real_frequencies .^ 2,
                   squared * azimuthal .^ 2,
                   vec(sum(2 .* abs2.(component_dr); dims=2)),
                   vec(sum(2 .* abs2.(component_dtheta); dims=2)),
                   -(squared * (real_frequencies .* azimuthal)),
                   phase_averaged_potential(quadrature, potential, profiles))
end

"""
    charge_densities(problem, s, bilinears)

Energy and angular-momentum content per unit tortoise coordinate at offset `s`,
integrated over the sphere. The `1/Δ` of the Boyer–Lindquist inverse metric
cancels against `dr/dr* = Δ/(r²+a²)`, so nothing is singular at the horizon.
"""
function charge_densities(problem::RadialProblem, s, bilinears::FieldBilinears)
    background, mass_coupling = problem.background, problem.mass_coupling
    grid = problem.quadrature.angular
    a, r = background.a, radius(background, s)
    horizon_function, sum_r2a2 = delta(background, s), r2a2(background, s)
    x, weights = grid.x, grid.weights
    sin2 = 1 .- x .^ 2

    tt_weight = (sum_r2a2^2 .- a^2 * horizon_function .* sin2) ./ sum_r2a2        # −ρ² g^{tt} Δ/(r²+a²)
    phiphi_weight = (horizon_function .- a^2 .* sin2) ./ (sin2 .* sum_r2a2)         #  ρ² g^{φφ} Δ/(r²+a²)
    rr_weight = horizon_function^2 / sum_r2a2                                       #  ρ² g^{rr} Δ/(r²+a²)
    thetatheta_weight = horizon_function / sum_r2a2                                 #  ρ² g^{θθ} Δ/(r²+a²)
    potential_weight = rho2.(Ref(background), s, x) .* (horizon_function / sum_r2a2)
    tphi_weight = -2a * r / sum_r2a2                                                #  ρ² g^{tφ} Δ/(r²+a²)

    energy = 2pi * sum(weights .* (tt_weight .* bilinears.dt_dt ./ 2 .+ phiphi_weight .* bilinears.dphi_dphi ./ 2 .+
                                   rr_weight .* bilinears.dr_dr ./ 2 .+ thetatheta_weight .* bilinears.dtheta_dtheta ./ 2 .+
                                   potential_weight .* mass_coupling^2 .* bilinears.potential))
    angular_momentum = 2pi * sum(weights .* (-tt_weight .* bilinears.dt_dphi .+ tphi_weight .* bilinears.dphi_dphi))
    return Charges(energy, angular_momentum)
end

"""
    simpson_weights(grid)

Composite Simpson weights on an arbitrary (possibly non-uniform) grid: the
quadratic through each consecutive triple of nodes is integrated over its two
intervals, and a trailing odd interval is integrated with the quadratic through
the last three nodes. Exact for cubics on a uniform grid, for quadratics otherwise.
"""
function simpson_weights(grid::AbstractVector)
    n = length(grid)
    weights = zeros(n)
    node = 1
    while node + 2 ≤ n
        left_width, right_width = grid[node + 1] - grid[node], grid[node + 2] - grid[node + 1]
        total = left_width + right_width
        weights[node] += total / 6 * (2 - right_width / left_width)
        weights[node + 1] += total^3 / (6left_width * right_width)
        weights[node + 2] += total / 6 * (2 - left_width / right_width)
        node += 2
    end
    if node + 1 == n   # one interval left: quadratic through the last three nodes over the last interval
        left_width, right_width = grid[n - 1] - grid[n - 2], grid[n] - grid[n - 1]
        total = left_width + right_width
        weights[n] += right_width * (2right_width + 3left_width) / (6total)
        weights[n - 1] += right_width * (right_width + 3left_width) / (6left_width)
        weights[n - 2] -= right_width^3 / (6left_width * total)
    end
    return weights
end

"""Physical node state (profiles and slopes multiplied by `scale`)."""
physical_state(layout::UnknownLayout, unknowns::AbstractVector, scale, node::Int) =
    (state = NodeState(unknowns[node_range(layout, node)]); NodeState(scale .* state.profiles, scale .* state.slopes))

"""Real component frequencies `ω_q` decoded from the unknown vector."""
function real_component_frequencies(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector)
    frequency_reals = unknowns[parameter_range(layout)][frequency_range(layout)]
    return real.(component_frequencies(problem.mode_set, _complexify(frequency_reals)))
end

"""
    cloud_charges(problem, layout, unknowns, scale)
    cloud_charges(solution)

Energy and angular momentum of the field between the numerical endpoints.
"""
function cloud_charges(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale)
    mesh = problem.mesh
    frequencies = real_component_frequencies(problem, layout, unknowns)
    weights = simpson_weights(mesh.rstar)
    total = Charges(zero(eltype(unknowns)), zero(eltype(unknowns)))
    for node in 1:layout.n_nodes
        state = physical_state(layout, unknowns, scale, node)
        bilinears = field_bilinears(problem, mesh.offsets[node], state.profiles, state.slopes, frequencies)
        total = total + weights[node] * charge_densities(problem, mesh.offsets[node], bilinears)
    end
    return total
end

cloud_charges(solution::RadialSolution) = cloud_charges(solution.problem, solution.layout, solution.unknowns, solution.scale)

"""
    component_charges(solution, component_index)

Free-field (quadratic) energy and angular momentum carried by one Fourier
component alone, with every other channel switched off. This is the
per-cloud energy used to compare clouds of different Fourier content; the
nonlinear potential energy is not attributed to any single component.
"""
function component_charges(solution::RadialSolution, component_index::Int)
    problem, layout, scale = solution.problem, solution.layout, solution.scale
    mesh, mode_set = problem.mesh, problem.mode_set
    channels = mode_set.component_channels[component_index]
    frequencies = real_component_frequencies(problem, layout, solution.unknowns)
    weights = simpson_weights(mesh.rstar)
    total = Charges(0.0, 0.0)
    for node in 1:layout.n_nodes
        state = physical_state(layout, solution.unknowns, scale, node)
        only_this = zeros(eltype(state.profiles), length(state.profiles))
        only_this_slopes = copy(only_this)
        only_this[channels] = state.profiles[channels]
        only_this_slopes[channels] = state.slopes[channels]
        bilinears = field_bilinears(problem, mesh.offsets[node], only_this, only_this_slopes, frequencies;
                                    potential=FreePotential())
        total = total + weights[node] * charge_densities(problem, mesh.offsets[node], bilinears)
    end
    return total
end

"""
    outward_flux(problem, layout, unknowns, scale, node)
    outward_flux(solution, node)

Outward energy and angular-momentum fluxes through the sphere at mesh `node`:
`F_E = 4π(r²+a²) Σ_q ω_q Σ_l Im[R*_{l,q} Q_{l,q}]` and the same with `m_q` for `F_J`.
"""
function outward_flux(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale, node::Int)
    mode_set = problem.mode_set
    state = physical_state(layout, unknowns, scale, node)
    frequencies = real_component_frequencies(problem, layout, unknowns)
    wronskian = imag.(conj.(state.profiles) .* state.slopes)
    prefactor = 4pi * r2a2(problem.background, problem.mesh.offsets[node])
    energy = prefactor * sum(frequencies[mode_set.channel_component[channel]] * wronskian[channel]
                             for channel in 1:n_channels(mode_set))
    angular_momentum = prefactor * sum(channel_m(mode_set, channel) * wronskian[channel]
                                       for channel in 1:n_channels(mode_set))
    return Charges(energy, angular_momentum)
end

outward_flux(solution::RadialSolution, node::Int) =
    outward_flux(solution.problem, solution.layout, solution.unknowns, solution.scale, node)

"""
    boundary_losses(solution)

Signed net losses `L = F(outer) − F(inner)`, so that `dE/dt = −L_E` and
`dJ/dt = −L_J` on a fixed background. Negative values mean the cloud gains.
"""
boundary_losses(solution::RadialSolution) =
    outward_flux(solution, solution.layout.n_nodes) - outward_flux(solution, 1)

"""
    charge_jacobian(solution)

The `2 × K` matrix `[∂E/∂A_i; ∂J/∂A_i]` including the response of the whole
coupled solution (profiles, frequencies, driven components) to each amplitude.
Computed exactly for the discretised problem via the implicit-function
theorem: `du/dA_i = −J⁻¹ ∂F/∂A_i` using the factorised Newton Jacobian, after
which the charges are differentiated along that direction with one dual pass.
"""
function charge_jacobian(solution::RadialSolution)
    problem, layout, unknowns, scale = solution.problem, solution.layout, solution.unknowns, solution.scale
    outer_node = layout.n_nodes
    outer_block = unknowns[node_range(layout, outer_node)]
    parameter_block = unknowns[parameter_range(layout)]
    # Only the outer boundary block depends explicitly on the amplitudes.
    residual_wrt_amplitudes = ForwardDiff.jacobian(solution.amplitudes) do amplitudes
        boundary_residual(problem, layout, :outer, outer_block, parameter_block, amplitudes, scale)
    end
    result = zeros(2, layout.n_modes)
    for mode_index in 1:layout.n_modes
        forcing = zeros(n_unknowns(layout))
        forcing[outer_rows(layout)] = -residual_wrt_amplitudes[:, mode_index]
        direction = solution.jacobian_factorization \ forcing
        result[:, mode_index] = ForwardDiff.derivative(0.0) do epsilon
            collect(cloud_charges(problem, layout, unknowns .+ epsilon .* direction, scale))
        end
    end
    return result
end
