"""
Phase-averaged conserved charges of the interacting field, the fluxes that
carry them through the numerical boundaries, and the response of the charges
to the prescribed amplitudes.

Energies are in units of `f_a² M` and angular momenta in `f_a² M²`; fluxes are
per unit `t/M`. With the Killing vectors `∂_t`, `∂_φ` of Kerr,

    E = −∫ T^t_t √−g d³x,    J = ∫ T^t_φ √−g d³x,

and the outward fluxes through a sphere are `F_E = −∮ T^r_t √−g`, `F_J = ∮ T^r_φ √−g`,
so that `dE/dt = −(F_E(out) − F_E(in))` and likewise for `J`.

Each independent phase `ξ_i` also carries a wave action `N_i = ∂⟨L⟩/∂ω_i`
(Whitham), the Noether charge of the phase shift `ξ_i → ξ_i + c` of the
averaged Lagrangian. Component `q` contributes its Klein–Gordon charge `n_q`
with weight `q_i`, so `N_i = Σ_q q_i n_q`, `J = Σ_i m_i N_i` exactly, and at
fixed background `dE = Σ_i ω_i dN_i`. The actions are adiabatic invariants: a
slowly changing black hole changes them only through boundary fluxes, which is
what fixes the metric work in the energy balance (see `evolution.jl`).
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
    potential_slope(potential, field_value)

`U'(Φ) = Φ − h(Φ)`, from the nonlinearity `h` that every potential defines.
"""
potential_slope(potential::AxionPotential, field_value) = field_value - nonlinearity(potential, field_value)

"""Allocated `FieldBilinears` over `n_polar` polar nodes."""
FieldBilinears{T}(n_polar::Int) where {T} = FieldBilinears{T}((zeros(T, n_polar) for _ in 1:6)...)

"""Value type of a possibly dual number type (the type itself for a plain number)."""
_value_type(::Type{T}) where {T<:Real} = T
_value_type(::Type{<:ForwardDiff.Dual{Tag, V}}) where {Tag, V} = V
_partial_count(::Type{<:Real}) = 0
_partial_count(::Type{<:ForwardDiff.Dual{Tag, V, N}}) where {Tag, V, N} = N

"""
    DensityWorkspace{T}(quadrature)

Scratch arrays for the field bilinears and charge densities at one radial node
(element type `T`, which may be a ForwardDiff dual number): the physical
profiles and slopes of all channels, the stacked component sums
`[Re; Im][component, polar node]` of the profiles and of the slopes, the
`FieldBilinears`, and for the potential average a value-type `FieldWorkspace`
with buffers for the value profiles, one direction's partial profiles and the
partials `[polar node, direction]`. One workspace serves one task at a time.
"""
struct DensityWorkspace{T<:Number, V<:Real, W<:FieldWorkspace{V}}
    profiles::Vector{Complex{T}}
    slopes::Vector{Complex{T}}
    component_values::Matrix{T}
    component_slopes::Matrix{T}
    bilinears::FieldBilinears{T}
    potential_field::W
    value_profiles::Vector{Complex{V}}
    partial_profiles::Vector{Complex{V}}
    potential_partials::Matrix{V}
end

function DensityWorkspace{T}(quadrature::FieldQuadrature) where {T}
    n_chan, n_polar = n_channels(quadrature.mode_set), length(quadrature.angular.x)
    n_stacked = 2n_components(quadrature.mode_set)
    V = _value_type(T)
    potential_field = FieldWorkspace{V}(quadrature)
    DensityWorkspace{T, V, typeof(potential_field)}(zeros(Complex{T}, n_chan), zeros(Complex{T}, n_chan),
                                                    zeros(T, n_stacked, n_polar), zeros(T, n_stacked, n_polar),
                                                    FieldBilinears{T}(n_polar), potential_field, zeros(Complex{V}, n_chan),
                                                    zeros(Complex{V}, n_chan), zeros(V, n_polar, _partial_count(T)))
end

"""
    averaged_potential!(averages, workspace, quadrature, potential, profiles)

`⟨U(Φ)⟩_ξ` at each polar node, like `phase_averaged_potential`, written into
`averages` with the scratch arrays of a `DensityWorkspace`. For a workspace of
ForwardDiff duals the field is reconstructed from the values only, and each
partial is the exact derivative `2 Re Σ_q dF_q conj(⟨e^{+iq·ξ} U'(Φ)⟩)` (with
`dF_q` the partial of the polar-summed profile of component `q`), which needs
no arithmetic on dual-number matrices.
"""
averaged_potential!(averages::AbstractVector, workspace::DensityWorkspace, quadrature::FieldQuadrature,
                    potential::AxionPotential, profiles::AbstractVector) =
    phase_averaged_potential!(averages, workspace.potential_field, quadrature, potential, profiles)

function averaged_potential!(averages::AbstractVector, workspace::DensityWorkspace{D}, quadrature::FieldQuadrature,
                             potential::AxionPotential, profiles::AbstractVector) where {Tag, V, N, D<:ForwardDiff.Dual{Tag, V, N}}
    plain, partials = workspace.potential_field, workspace.potential_partials
    n_comp = n_components(quadrature.mode_set)
    workspace.value_profiles .= _value.(profiles)
    field = reconstruct_field!(plain, quadrature, workspace.value_profiles)
    plain.grid_values .= potential_slope.(Ref(potential), field)
    slope_average = project_phases!(plain.projected, quadrature, plain.grid_values)   # [Re; Im] ⟨e^{+iq·ξ} U′(Φ)⟩
    for slot in 1:N
        workspace.partial_profiles .= _partial.(profiles, slot)
        # The field is built, so the component buffer can take the partial profiles' component sums.
        partial_components = component_profiles!(plain.component_values, quadrature, workspace.partial_profiles)
        for node in axes(partials, 1)
            total = zero(V)
            for component_index in 1:n_comp
                total += partial_components[component_index, node] * slope_average[component_index, node] +
                         partial_components[n_comp + component_index, node] * slope_average[n_comp + component_index, node]
            end
            partials[node, slot] = 2total
        end
    end
    plain.grid_values .= potential_energy.(Ref(potential), field)
    for node in eachindex(averages)
        value = sum(view(plain.grid_values, :, node)) / n_points(quadrature.phase)
        averages[node] = D(value, ForwardDiff.Partials(ntuple(slot -> partials[node, slot], Val(N))))
    end
    return averages
end

"""
    averaged_potential(quadrature, potential, profiles)

Allocating `averaged_potential!`: `⟨U(Φ)⟩_ξ` at each polar node, with exact partials for dual-number profiles.
"""
function averaged_potential(quadrature::FieldQuadrature, potential::AxionPotential, profiles::AbstractVector)
    T = field_eltype(quadrature, profiles)
    return averaged_potential!(zeros(T, length(quadrature.angular.x)), DensityWorkspace{T}(quadrature), quadrature,
                               potential, profiles)
end

"""
    field_bilinears!(bilinears, workspace, problem, s, profiles, slopes, real_frequencies; potential=problem.potential)

Bilinears at offset `s` from the *physical* profiles and tortoise slopes of all
channels and the real component frequencies `ω_q`, written into `bilinears`
with the scratch arrays of a `DensityWorkspace`.
"""
function field_bilinears!(bilinears::FieldBilinears, workspace::DensityWorkspace, problem::RadialProblem, s,
                          profiles::AbstractVector, slopes::AbstractVector, real_frequencies::AbstractVector;
                          potential::AxionPotential=problem.potential)
    quadrature, mode_set = problem.quadrature, problem.mode_set
    harmonics_dtheta = quadrature.angular.harmonics_dtheta
    n_comp = n_components(mode_set)
    values = component_profiles!(workspace.component_values, quadrature, profiles)
    radial = component_profiles!(workspace.component_slopes, quadrature, slopes)   # ∂F_q/∂r* (stacked)
    radial_factor = dr_drstar(problem.background, s)
    T = eltype(bilinears.dt_dt)
    for node in axes(values, 2)
        dt_dt = dphi_dphi = dr_dr = dtheta_dtheta = dt_dphi = zero(T)
        for (component_index, component) in enumerate(mode_set.components)
            squared = 2 * (values[component_index, node]^2 + values[n_comp + component_index, node]^2)   # 2|F_q|²
            frequency, m = real_frequencies[component_index], component.m
            dt_dt += squared * frequency^2
            dphi_dphi += squared * m^2
            dt_dphi -= squared * (frequency * m)
            dr_dr += 2 * ((radial[component_index, node] / radial_factor)^2 +
                          (radial[n_comp + component_index, node] / radial_factor)^2)
            real_dtheta = imag_dtheta = zero(T)
            for channel in mode_set.component_channels[component_index]
                real_dtheta += harmonics_dtheta[node, channel] * real(profiles[channel])
                imag_dtheta += harmonics_dtheta[node, channel] * imag(profiles[channel])
            end
            dtheta_dtheta += 2 * (real_dtheta^2 + imag_dtheta^2)
        end
        bilinears.dt_dt[node], bilinears.dphi_dphi[node], bilinears.dr_dr[node] = dt_dt, dphi_dphi, dr_dr
        bilinears.dtheta_dtheta[node], bilinears.dt_dphi[node] = dtheta_dtheta, dt_dphi
    end
    averaged_potential!(bilinears.potential, workspace, quadrature, potential, profiles)
    return bilinears
end

"""
    field_bilinears(problem, s, profiles, slopes, real_frequencies; potential=problem.potential)

Allocating `field_bilinears!`.
"""
function field_bilinears(problem::RadialProblem, s, profiles::AbstractVector, slopes::AbstractVector,
                         real_frequencies::AbstractVector; potential::AxionPotential=problem.potential)
    T = promote_type(field_eltype(problem.quadrature, profiles, s, problem.background.a), real(eltype(slopes)),
                     eltype(real_frequencies))
    return field_bilinears!(FieldBilinears{T}(length(problem.quadrature.angular.x)), DensityWorkspace{T}(problem.quadrature),
                            problem, s, profiles, slopes, real_frequencies; potential)
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
    rr_weight = horizon_function^2 / sum_r2a2                  #  ρ² g^{rr} Δ/(r²+a²)
    thetatheta_weight = horizon_function / sum_r2a2            #  ρ² g^{θθ} Δ/(r²+a²)
    tphi_weight = -2a * r / sum_r2a2                           #  ρ² g^{tφ} Δ/(r²+a²)
    energy = angular_momentum = zero(promote_type(eltype(bilinears.dt_dt), typeof(rr_weight), typeof(mass_coupling)))
    for (node, x) in enumerate(grid.x)
        sin2 = 1 - x^2
        tt_weight = (sum_r2a2^2 - a^2 * horizon_function * sin2) / sum_r2a2       # −ρ² g^{tt} Δ/(r²+a²)
        phiphi_weight = (horizon_function - a^2 * sin2) / (sin2 * sum_r2a2)        #  ρ² g^{φφ} Δ/(r²+a²)
        potential_weight = rho2(background, s, x) * (horizon_function / sum_r2a2)
        energy += grid.weights[node] * (tt_weight * bilinears.dt_dt[node] / 2 + phiphi_weight * bilinears.dphi_dphi[node] / 2 +
                                        rr_weight * bilinears.dr_dr[node] / 2 + thetatheta_weight * bilinears.dtheta_dtheta[node] / 2 +
                                        potential_weight * mass_coupling^2 * bilinears.potential[node])
        angular_momentum += grid.weights[node] * (-tt_weight * bilinears.dt_dphi[node] + tphi_weight * bilinears.dphi_dphi[node])
    end
    return Charges(2pi * energy, 2pi * angular_momentum)
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
    integrated_densities(problem, layout, unknowns, scale; charges=true, actions=true)

Simpson integrals over the mesh of the per-node densities, as the vector
`[E, J]` (with `charges`) followed by the Klein–Gordon charge `n_q` of every
component (with `actions`). Nodes are independent and are evaluated on the
available threads into one column each; the sum then runs in node order, so the
result does not depend on the thread count.
"""
function integrated_densities(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale;
                              charges::Bool=true, actions::Bool=true)
    mesh = problem.mesh
    frequencies = real_component_frequencies(problem, layout, unknowns)
    weights = simpson_weights(ForwardDiff.value.(mesh.rstar))
    row_count = (charges ? 2 : 0) + (actions ? n_components(problem.mode_set) : 0)
    element_type = promote_type(eltype(unknowns), typeof(scale), typeof(problem.mass_coupling), eltype(mesh.offsets),
                                typeof(problem.background.a))
    densities = Matrix{element_type}(undef, row_count, layout.n_nodes)
    _foreach_in_chunks(1:layout.n_nodes, () -> DensityWorkspace{element_type}(problem.quadrature)) do workspace, node
        _node_densities!(view(densities, :, node), workspace, problem, layout, unknowns, scale, node, frequencies;
                         charges, actions)
    end
    totals = zeros(element_type, row_count)
    for node in 1:layout.n_nodes
        totals .+= weights[node] .* view(densities, :, node)
    end
    return totals
end

"""
Densities of `integrated_densities` at mesh `node` (`[E, J]` with `charges`, then `n_q` with `actions`), written
into `column` with the scratch arrays of a `DensityWorkspace`.
"""
function _node_densities!(column::AbstractVector, workspace::DensityWorkspace, problem::RadialProblem, layout::UnknownLayout,
                          unknowns::AbstractVector, scale, node::Int, frequencies::AbstractVector; charges::Bool, actions::Bool)
    block, offset = view(unknowns, node_range(layout, node)), problem.mesh.offsets[node]
    for channel in eachindex(workspace.profiles)
        workspace.profiles[channel] = scale * complex(block[4channel - 3], block[4channel - 2])
        workspace.slopes[channel] = scale * complex(block[4channel - 1], block[4channel])
    end
    row = 0
    if charges
        bilinears = field_bilinears!(workspace.bilinears, workspace, problem, offset, workspace.profiles, workspace.slopes,
                                     frequencies)
        node_charges = charge_densities(problem, offset, bilinears)
        column[1], column[2] = node_charges.energy, node_charges.angular_momentum
        row = 2
    end
    actions && component_action_densities!(view(column, row + 1:length(column)), workspace.component_values, problem, offset,
                                            workspace.profiles, frequencies)
    return column
end

"""
    cloud_charges(problem, layout, unknowns, scale)
    cloud_charges(solution)

Energy and angular momentum of the field between the numerical endpoints.
"""
function cloud_charges(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale)
    totals = integrated_densities(problem, layout, unknowns, scale; actions=false)
    return Charges(totals[1], totals[2])
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
    component_action_densities!(densities, component_values, problem, s, profiles, real_frequencies)

Klein–Gordon charge of every Fourier component per unit tortoise coordinate at
offset `s`, `dn_q/dr* = 2π Σ_x w ρ² (Δ/(r²+a²)) 2|F_q|² (−g^{tt} ω_q + g^{tφ} m_q)`,
from the physical profiles and the real component frequencies, written into
`densities`; `component_values` is scratch for the stacked component profiles.
"""
function component_action_densities!(densities::AbstractVector, component_values::AbstractMatrix, problem::RadialProblem, s,
                                     profiles::AbstractVector, real_frequencies::AbstractVector)
    background, quadrature, mode_set = problem.background, problem.quadrature, problem.mode_set
    grid = quadrature.angular
    n_comp = n_components(mode_set)
    a, r = background.a, radius(background, s)
    horizon_function, sum_r2a2 = delta(background, s), r2a2(background, s)
    tphi_weight = -2a * r / sum_r2a2                                         #  ρ² g^{tφ} Δ/(r²+a²)
    values = component_profiles!(component_values, quadrature, profiles)
    for (component_index, component) in enumerate(mode_set.components)
        total = zero(eltype(densities))
        for (node, x) in enumerate(grid.x)
            tt_weight = (sum_r2a2^2 - a^2 * horizon_function * (1 - x^2)) / sum_r2a2   # −ρ² g^{tt} Δ/(r²+a²)
            squared = 2 * (values[component_index, node]^2 + values[n_comp + component_index, node]^2)   # 2|F_q|²
            total += grid.weights[node] * squared * (tt_weight * real_frequencies[component_index] + tphi_weight * component.m)
        end
        densities[component_index] = 2pi * total
    end
    return densities
end

"""
    component_action_densities(problem, s, profiles, real_frequencies)

Allocating `component_action_densities!`.
"""
function component_action_densities(problem::RadialProblem, s, profiles::AbstractVector, real_frequencies::AbstractVector)
    T = promote_type(field_eltype(problem.quadrature, profiles, s, problem.background.a), eltype(real_frequencies))
    n_stacked, n_polar = 2n_components(problem.mode_set), length(problem.quadrature.angular.x)
    return component_action_densities!(zeros(T, n_components(problem.mode_set)), zeros(T, n_stacked, n_polar), problem, s,
                                       profiles, real_frequencies)
end

"""Matrix `[mode, component]` of the integer vectors `q`: `N_i = Σ_q q_i n_q`."""
phase_weights(mode_set::ModeSet) = [component.q[mode] for mode in 1:n_modes(mode_set), component in mode_set.components]

"""
    phase_actions(problem, layout, unknowns, scale)
    phase_actions(solution)

Wave actions `N_i = Σ_q q_i n_q` of the independent phases between the
numerical endpoints, in units of `f_a² M²` (like `J`), as a vector over modes.
"""
function phase_actions(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, scale)
    return phase_weights(problem.mode_set) * integrated_densities(problem, layout, unknowns, scale; charges=false)
end

phase_actions(solution::RadialSolution) = phase_actions(solution.problem, solution.layout, solution.unknowns, solution.scale)

"""
    action_flux(solution, node)

Outward action fluxes of the independent phases through the sphere at mesh
`node`, `F_{N_i} = Σ_q q_i 4π(r²+a²) Σ_l Im[R*_{l,q} Q_{l,q}]`, so that a
component's energy and angular-momentum fluxes are `ω_q` and `m_q` times its
action flux (`outward_flux`).
"""
function action_flux(solution::RadialSolution, node::Int)
    problem, layout = solution.problem, solution.layout
    mode_set = problem.mode_set
    state = physical_state(layout, solution.unknowns, solution.scale, node)
    wronskian = imag.(conj.(state.profiles) .* state.slopes)
    prefactor = 4pi * r2a2(problem.background, problem.mesh.offsets[node])
    component_fluxes = [prefactor * sum(wronskian[channels]) for channels in mode_set.component_channels]
    return phase_weights(mode_set) * component_fluxes
end

"""
    action_losses(solution)

Signed net action losses `F_N(outer) − F_N(inner)` of the independent phases,
so that `dN_i/dt = −L_{N_i}` on any slowly varying background.
"""
action_losses(solution::RadialSolution) = action_flux(solution, solution.layout.n_nodes) - action_flux(solution, 1)

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
    amplitude_response(solution)

Directions `du/dA_i` (columns) of the whole discretised solution with respect
to each amplitude, by the implicit-function theorem `du/dA_i = −J⁻¹ ∂F/∂A_i`
with the factorised Newton Jacobian. Only the outer boundary block depends
explicitly on the amplitudes.
"""
function amplitude_response(solution::RadialSolution)
    problem, layout, unknowns, scale = solution.problem, solution.layout, solution.unknowns, solution.scale
    outer_block = unknowns[node_range(layout, layout.n_nodes)]
    parameter_block = unknowns[parameter_range(layout)]
    residual_wrt_amplitudes = ForwardDiff.jacobian(solution.amplitudes) do amplitudes
        boundary_residual(problem, layout, :outer, outer_block, parameter_block, amplitudes, scale)
    end
    forcing = zeros(n_unknowns(layout), layout.n_modes)
    forcing[outer_rows(layout), :] = -residual_wrt_amplitudes
    return solution.jacobian_factorization \ forcing   # one back-substitution with K right-hand sides
end

"""
    observable_jacobian(solution, observable; directions=amplitude_response(solution))

Matrix `[∂O_k/∂A_i]` of a vector-valued `observable(problem, layout, unknowns, scale)`
including the response of the coupled solution, with one dual pass of the
observable along each amplitude direction.
"""
function observable_jacobian(solution::RadialSolution, observable; directions=amplitude_response(solution))
    problem, layout, unknowns, scale = solution.problem, solution.layout, solution.unknowns, solution.scale
    columns = map(eachcol(directions)) do direction
        ForwardDiff.derivative(epsilon -> observable(problem, layout, unknowns .+ epsilon .* direction, scale), 0.0)
    end
    return reduce(hcat, columns)
end

"""
    charge_jacobian(solution)

The `2 × K` matrix `[∂E/∂A_i; ∂J/∂A_i]` including the response of the whole
coupled solution (profiles, frequencies, driven components) to each amplitude.
Computed exactly for the discretised problem via the implicit-function
theorem (`amplitude_response`), after which the charges are differentiated
along each direction with one dual pass.
"""
charge_jacobian(solution::RadialSolution; directions=amplitude_response(solution)) =
    observable_jacobian(solution, (arguments...) -> collect(cloud_charges(arguments...)); directions)

"""
    action_jacobian(solution)

The `K × K` matrix `[∂N_j/∂A_i]` of the phase actions, like `charge_jacobian`.
"""
action_jacobian(solution::RadialSolution; directions=amplitude_response(solution)) =
    observable_jacobian(solution, phase_actions; directions)
