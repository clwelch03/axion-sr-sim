"""
Reconstruction of the real field on a `(x, ξ)` quadrature grid and its
projection back onto `(l, q)` channels. This is the one place that knows how
to turn radial profiles into a physical field, and it serves both the
nonlinear source `N_{l,q}` and the phase-averaged potential energy.
"""

"""
    PhaseGrid(n_per_dim, mode_set)

Uniform tensor grid on the torus `[0, 2π)^K` of independent phases, with the
table `phase_factors[point, component] = exp(-i q · ξ_point)`. The trapezoidal
rule on a periodic grid is spectrally accurate, and `n_per_dim` must exceed
the highest harmonic that the nonlinearity can alias back onto a retained
component.

The same table is also stored in real form, so that synthesis and projection
run as real matrix products on stacked real and imaginary parts (for `K′`
components, the first `K′` rows of a stacked matrix hold real parts and the next
`K′` imaginary parts):

* `synthesis[point, :] = [2cos(q·ξ_point), 2sin(q·ξ_point)]`, so that
  `Φ = synthesis * [Re F; Im F]`;
* `analysis[:, point] = [cos(q·ξ_point); sin(q·ξ_point)] / n_points`, so that
  `analysis * f = [Re; Im] ⟨e^{+i q·ξ} f⟩_ξ`.

For the source Jacobian it also stores the averaging table
`exp(+i k_slot · ξ_point) / n_points` in the same stacked form,
`jacobian_table = [cos(k_slot·ξ_point); sin(k_slot·ξ_point)] / n_points`, for the
distinct wavevectors `k = q ± q′` over all component pairs up to sign: `h′(Φ)` is
real, so `ĥ′(−k) = conj(ĥ′(k))` and each `±k` pair needs one row.
`sum_slot[q, q′]`, `difference_slot[q, q′]` point into the table and
`sum_conjugated`, `difference_conjugated` say whether the stored row is for
`−k`, so that every phase average the Jacobian needs comes from one matrix
product (see `nonlinear_source_jacobian`).
"""
struct PhaseGrid{T<:Real}
    n_per_dim::Int
    phase_factors::Matrix{Complex{T}}
    synthesis::Matrix{T}         # [point, stacked component]
    analysis::Matrix{T}          # [stacked component, point]
    jacobian_table::Matrix{T}    # [stacked slot, point]
    sum_slot::Matrix{Int}
    difference_slot::Matrix{Int}
    sum_conjugated::Matrix{Bool}
    difference_conjugated::Matrix{Bool}
end

function PhaseGrid(n_per_dim::Int, mode_set::ModeSet)
    one_axis = 2pi .* (0:n_per_dim-1) ./ n_per_dim
    torus_points = vec(collect(Iterators.product(ntuple(_ -> one_axis, n_modes(mode_set))...)))
    phase_factors = [cis(-sum(component.q .* collect(point)))
                     for point in torus_points, component in mode_set.components]

    # Distinct wavevectors up to sign; a slot stores the representative with the larger lexicographic order.
    wavevectors = Vector{Vector{Int}}()
    function slot_of(wavevector)
        representative = max(wavevector, -wavevector)
        slot = findfirst(==(representative), wavevectors)
        if slot === nothing
            push!(wavevectors, representative)
            slot = length(wavevectors)
        end
        return slot
    end
    components = mode_set.components
    sum_slot = [slot_of(target.q .+ source.q) for target in components, source in components]
    difference_slot = [slot_of(target.q .- source.q) for target in components, source in components]
    sum_conjugated = [target.q .+ source.q != max(target.q .+ source.q, -(target.q .+ source.q))
                      for target in components, source in components]
    difference_conjugated = [target.q .- source.q != max(target.q .- source.q, -(target.q .- source.q))
                             for target in components, source in components]
    averaging_table = [cis(sum(wavevector .* collect(point))) / length(torus_points)
                       for wavevector in wavevectors, point in torus_points]
    # cos(q·ξ) = Re exp(-iq·ξ), sin(q·ξ) = -Im exp(-iq·ξ)
    synthesis = hcat(2 .* real.(phase_factors), -2 .* imag.(phase_factors))
    analysis = vcat(transpose(real.(phase_factors)), -transpose(imag.(phase_factors))) ./ length(torus_points)
    PhaseGrid(n_per_dim, phase_factors, synthesis, analysis, vcat(real.(averaging_table), imag.(averaging_table)),
              sum_slot, difference_slot, sum_conjugated, difference_conjugated)
end

n_points(grid::PhaseGrid) = size(grid.phase_factors, 1)
n_slots(grid::PhaseGrid) = size(grid.jacobian_table, 1) ÷ 2

"""
    FieldQuadrature(mode_set; n_polar=12, n_phase=16)

Everything needed to evaluate phase- and polar-averaged nonlinear quantities:
the polar Gauss–Legendre grid, the phase torus grid, and the mode bookkeeping.
"""
struct FieldQuadrature{T<:Real}
    mode_set::ModeSet
    angular::AngularGrid{T}
    phase::PhaseGrid{T}
end

FieldQuadrature(mode_set::ModeSet; n_polar::Int=12, n_phase::Int=16) =
    FieldQuadrature(mode_set, AngularGrid(n_polar, mode_set), PhaseGrid(n_phase, mode_set))

"""
    FieldWorkspace{T}(quadrature)

Scratch arrays for the in-place field kernels at one radial node (element type
`T`, which may be a dual number): the stacked component profiles
`component_values = [Re F; Im F][component, polar node]`, the field
`Φ[point, polar node]`, a grid function of the field (`h(Φ)`, `h′(Φ)` or
`U(Φ)`), its stacked phase projection `[Re; Im] ⟨e^{+i q·ξ} h⟩[component, polar node]`,
the stacked Fourier coefficients `[Re; Im] ĥ′(k)[slot, polar node]`, and the polar
weights `w ρ²`. One workspace serves one task at a time.
"""
struct FieldWorkspace{T<:Number}
    component_values::Matrix{T}
    field::Matrix{T}
    grid_values::Matrix{T}
    projected::Matrix{T}
    fourier_coefficients::Matrix{T}
    polar_weights::Vector{T}
end

function FieldWorkspace{T}(quadrature::FieldQuadrature) where {T}
    n_polar, n_stacked = length(quadrature.angular.x), 2n_components(quadrature.mode_set)
    n_grid = n_points(quadrature.phase)
    FieldWorkspace{T}(zeros(T, n_stacked, n_polar), zeros(T, n_grid, n_polar), zeros(T, n_grid, n_polar),
                      zeros(T, n_stacked, n_polar), zeros(T, 2n_slots(quadrature.phase), n_polar), zeros(T, n_polar))
end

"""
    component_profiles!(component_values, quadrature, profiles)

Stacked component profiles `[Re F_q; Im F_q][component, polar node]` with
`F_q(x_node) = Σ_l R_{l,q} Ỹ_{l,|m_q|}(x_node)`, written into `component_values`.
"""
function component_profiles!(component_values::AbstractMatrix, quadrature::FieldQuadrature, profiles::AbstractVector)
    mode_set, harmonics = quadrature.mode_set, quadrature.angular.harmonics
    n_comp = n_components(mode_set)
    for (component_index, channels) in enumerate(mode_set.component_channels), node in axes(harmonics, 1)
        real_part = imag_part = zero(eltype(component_values))
        for channel in channels
            real_part += harmonics[node, channel] * real(profiles[channel])
            imag_part += harmonics[node, channel] * imag(profiles[channel])
        end
        component_values[component_index, node] = real_part
        component_values[n_comp + component_index, node] = imag_part
    end
    return component_values
end

"""
    component_profiles(quadrature, profiles)

Sum each component's channels over their polar harmonics, returning
`F[node, component] = F_q(x_node) = Σ_l R_{l,q} Ỹ_{l,|m_q|}(x_node)`.
"""
function component_profiles(quadrature::FieldQuadrature, profiles::AbstractVector)
    n_comp = n_components(quadrature.mode_set)
    T = promote_type(real(eltype(profiles)), eltype(quadrature.angular.harmonics))
    stacked = component_profiles!(zeros(T, 2n_comp, length(quadrature.angular.x)), quadrature, profiles)
    return permutedims(complex.(stacked[1:n_comp, :], stacked[n_comp+1:end, :]))
end

"""
    reconstruct_field!(workspace, quadrature, profiles)

Fill `workspace.field` with `Φ(x_node, ξ_point) = Σ_q F_q(x_node) e^{-i q·ξ_point} + c.c.`
(indexed `[point, node]`) and return it; `workspace.component_values` holds the
stacked component profiles afterwards.
"""
function reconstruct_field!(workspace::FieldWorkspace, quadrature::FieldQuadrature, profiles::AbstractVector)
    component_profiles!(workspace.component_values, quadrature, profiles)
    return mul!(workspace.field, quadrature.phase.synthesis, workspace.component_values)
end

"""Element type of field workspaces for these profiles and extra scalars (dual numbers propagate)."""
field_eltype(quadrature::FieldQuadrature, profiles::AbstractVector, scalars...) =
    promote_type(real(eltype(profiles)), eltype(quadrature.angular.x), map(typeof, scalars)...)

"""
    reconstruct_field(quadrature, profiles)

The real field `Φ(x_node, ξ_point) = Σ_q F_q(x_node) e^{-i q·ξ_point} + c.c.` as a
matrix indexed `[point, node]` (phase point, polar node).
"""
function reconstruct_field(quadrature::FieldQuadrature, profiles::AbstractVector)
    workspace = FieldWorkspace{field_eltype(quadrature, profiles)}(quadrature)
    return reconstruct_field!(workspace, quadrature, profiles)
end

"""
    project_phases!(projected, quadrature, grid_function)

Stacked phase average `[Re; Im] ⟨e^{+i q·ξ} f⟩_ξ[component, node]` of a grid
function `f[point, node]`, written into `projected`.
"""
project_phases!(projected::AbstractMatrix, quadrature::FieldQuadrature, grid_function::AbstractMatrix) =
    mul!(projected, quadrature.phase.analysis, grid_function)

"""
    project_phases(quadrature, grid_function)

Phase average `⟨e^{+i q·ξ} f⟩_ξ` of a grid function `f[point, node]`, returning
a matrix indexed `[component, node]`.
"""
function project_phases(quadrature::FieldQuadrature, grid_function::AbstractMatrix)
    n_comp = n_components(quadrature.mode_set)
    stacked = quadrature.phase.analysis * grid_function
    return complex.(stacked[1:n_comp, :], stacked[n_comp+1:end, :])
end

"""Fill `workspace.polar_weights` with the Gauss weights times `ρ²(r, x)` at offset `s`."""
function _polar_weights!(workspace::FieldWorkspace, quadrature::FieldQuadrature, background::KerrBackground, s)
    grid = quadrature.angular
    workspace.polar_weights .= grid.weights .* rho2.(Ref(background), s, grid.x)
    return workspace.polar_weights
end

"""
    nonlinear_source!(source, workspace, quadrature, potential, background, α, s, profiles)

In-place `nonlinear_source`: writes `N_{l,q}` of every channel into `source`
(complex vector) using the scratch arrays of `workspace`.
"""
function nonlinear_source!(source::AbstractVector, workspace::FieldWorkspace, quadrature::FieldQuadrature,
                           potential::AxionPotential, background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    mode_set, harmonics = quadrature.mode_set, quadrature.angular.harmonics
    n_comp = n_components(mode_set)
    field = reconstruct_field!(workspace, quadrature, profiles)
    workspace.grid_values .= nonlinearity.(Ref(potential), field)
    projected = project_phases!(workspace.projected, quadrature, workspace.grid_values)
    weights = _polar_weights!(workspace, quadrature, background, s)
    for channel in eachindex(source)
        component_index = mode_set.channel_component[channel]
        real_part = imag_part = zero(eltype(weights))
        for node in eachindex(weights)
            weighted_harmonic = weights[node] * harmonics[node, channel]
            real_part += weighted_harmonic * projected[component_index, node]
            imag_part += weighted_harmonic * projected[n_comp + component_index, node]
        end
        source[channel] = complex(real_part, imag_part) * mass_coupling^2
    end
    return source
end

nonlinear_source!(source::AbstractVector, ::FieldWorkspace, ::FieldQuadrature, ::FreePotential, ::KerrBackground,
                  mass_coupling, s, profiles::AbstractVector) = fill!(source, zero(eltype(source)))

"""
    nonlinear_source(quadrature, potential, background, α, s, profiles)

The projected nonlinear source for every channel at one radial node,
`N_{l,q} = α² ∫ dx Ỹ_{l,|m_q|}(x) ρ²(r, x) ⟨e^{+i q·ξ} h(Φ)⟩_ξ`,
evaluated on the full reconstructed field without any expansion in the
self-interaction strength.
"""
function nonlinear_source(quadrature::FieldQuadrature, potential::AxionPotential,
                          background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    T = field_eltype(quadrature, profiles, mass_coupling, s, background.a)
    source = Vector{Complex{T}}(undef, length(profiles))
    return nonlinear_source!(source, FieldWorkspace{T}(quadrature), quadrature, potential, background, mass_coupling, s, profiles)
end

nonlinear_source(quadrature::FieldQuadrature, ::FreePotential, ::KerrBackground, mass_coupling, s, profiles::AbstractVector) =
    zeros(promote_type(eltype(profiles), typeof(mass_coupling)), length(profiles))

"""
    NonlinearSourceJacobian

Derivative of the nonlinear source at one radial node with respect to the
real and imaginary parts of every channel's profile: `wrt_real[c, c′] =
∂N_c/∂(Re R_{c′})` and `wrt_imag[c, c′] = ∂N_c/∂(Im R_{c′})`, both complex.
The source is real-differentiable but not complex-linear, so the two matrices
are independent.
"""
struct NonlinearSourceJacobian{T<:Real}
    wrt_real::Matrix{Complex{T}}
    wrt_imag::Matrix{Complex{T}}
end

NonlinearSourceJacobian{T}(n_chan::Int) where {T} =
    NonlinearSourceJacobian{T}(zeros(Complex{T}, n_chan, n_chan), zeros(Complex{T}, n_chan, n_chan))

"""
    nonlinear_source_jacobian!(jacobian, workspace, quadrature, potential, background, α, s, profiles)

In-place `nonlinear_source_jacobian`: overwrites `jacobian.wrt_real` and
`jacobian.wrt_imag` using the scratch arrays of `workspace`.
"""
function nonlinear_source_jacobian!(jacobian::NonlinearSourceJacobian, workspace::FieldWorkspace,
                                    quadrature::FieldQuadrature, potential::AxionPotential,
                                    background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    mode_set, grid, phase = quadrature.mode_set, quadrature.angular, quadrature.phase
    wrt_real, wrt_imag = jacobian.wrt_real, jacobian.wrt_imag
    fill!(wrt_real, 0)
    fill!(wrt_imag, 0)
    field = reconstruct_field!(workspace, quadrature, profiles)
    workspace.grid_values .= nonlinearity_derivative.(Ref(potential), field)   # h′(Φ)[point, node]
    fourier = mul!(workspace.fourier_coefficients, phase.jacobian_table, workspace.grid_values)  # [Re; Im] ĥ′(k)[slot, node]
    polar_weights = _polar_weights!(workspace, quadrature, background, s)
    polar_weights .*= mass_coupling^2
    slot_count = n_slots(phase)
    owner = mode_set.channel_component
    n_chan = size(wrt_real, 1)
    for node in eachindex(polar_weights), source in 1:n_chan, target in 1:n_chan
        target_component, source_component = owner[target], owner[source]
        sum_slot = phase.sum_slot[target_component, source_component]
        difference_slot = phase.difference_slot[target_component, source_component]
        plus = complex(fourier[sum_slot, node], fourier[slot_count + sum_slot, node])
        minus = complex(fourier[difference_slot, node], fourier[slot_count + difference_slot, node])
        phase.sum_conjugated[target_component, source_component] && (plus = conj(plus))
        phase.difference_conjugated[target_component, source_component] && (minus = conj(minus))
        prefactor = polar_weights[node] * grid.harmonics[node, target] * grid.harmonics[node, source]
        wrt_real[target, source] += prefactor * (plus + minus)
        wrt_imag[target, source] += prefactor * (-im) * (plus - minus)
    end
    return jacobian
end

function nonlinear_source_jacobian!(jacobian::NonlinearSourceJacobian, ::FieldWorkspace, ::FieldQuadrature, ::FreePotential,
                                    ::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    fill!(jacobian.wrt_real, 0)
    fill!(jacobian.wrt_imag, 0)
    return jacobian
end

"""
    nonlinear_source_jacobian(quadrature, potential, background, α, s, profiles)

Jacobian of `nonlinear_source` with respect to the profiles. With
`Φ = Σ_q 2 Re[F_q e^{-i q·ξ}]` one has `∂Φ/∂Re R_c = 2 cos(q_c·ξ) Ỹ_c` and
`∂Φ/∂Im R_c = 2 sin(q_c·ξ) Ỹ_c`, so

    ∂N_c/∂Re R_{c′} = α² Σ_x w ρ² Ỹ_c Ỹ_{c′} ⟨e^{+i q_c·ξ} h′(Φ) 2cos(q_{c′}·ξ)⟩_ξ,

and likewise with `2 sin` for the imaginary part. Writing the cosine and sine
as exponentials, both averages are combinations of the Fourier coefficients
`ĥ′(k) = ⟨e^{+i k·ξ} h′(Φ)⟩_ξ` at `k = q_c ± q_{c′}`:

    ⟨e^{i q·ξ} h′ 2cos(q′·ξ)⟩ = ĥ′(q + q′) + ĥ′(q − q′),
    ⟨e^{i q·ξ} h′ 2sin(q′·ξ)⟩ = −i [ĥ′(q + q′) − ĥ′(q − q′)].

The coefficients at all polar nodes come from a single real product of the
phase grid's stacked averaging table with `h′(Φ)[point, node]`.
"""
function nonlinear_source_jacobian(quadrature::FieldQuadrature, potential::AxionPotential,
                                   background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    T = field_eltype(quadrature, profiles, mass_coupling, s, background.a)
    return nonlinear_source_jacobian!(NonlinearSourceJacobian{T}(length(profiles)), FieldWorkspace{T}(quadrature),
                                      quadrature, potential, background, mass_coupling, s, profiles)
end

"""
    phase_averaged_potential!(averages, workspace, quadrature, potential, profiles)

In-place `phase_averaged_potential`: writes `⟨U(Φ)⟩_ξ` at each polar node into `averages`.
"""
function phase_averaged_potential!(averages::AbstractVector, workspace::FieldWorkspace, quadrature::FieldQuadrature,
                                   potential::AxionPotential, profiles::AbstractVector)
    field = reconstruct_field!(workspace, quadrature, profiles)
    workspace.grid_values .= potential_energy.(Ref(potential), field)
    sum!(transpose(averages), workspace.grid_values)
    averages ./= n_points(quadrature.phase)
    return averages
end

"""
    phase_averaged_potential(quadrature, potential, profiles)

`⟨U(Φ)⟩_ξ` at each polar node for one radial node, as a vector over nodes.
"""
function phase_averaged_potential(quadrature::FieldQuadrature, potential::AxionPotential, profiles::AbstractVector)
    T = field_eltype(quadrature, profiles)
    return phase_averaged_potential!(zeros(T, length(quadrature.angular.x)), FieldWorkspace{T}(quadrature),
                                     quadrature, potential, profiles)
end
