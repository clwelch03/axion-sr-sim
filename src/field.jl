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

For the source Jacobian it also stores the averaging table
`exp(+i k_slot · ξ_point) / n_points`, split into real matrices
`jacobian_cosines` and `jacobian_sines` (so the product with the real `h′(Φ)`
runs as BLAS), for the distinct wavevectors `k = q ± q′` over all component
pairs up to sign: `h′(Φ)` is real, so `ĥ′(−k) = conj(ĥ′(k))` and each `±k` pair
needs one row. `sum_slot[q, q′]`, `difference_slot[q, q′]` point into the table
and `sum_conjugated`, `difference_conjugated` say whether the stored row is for
`−k`, so that every phase average the Jacobian needs comes from one matrix
product (see `nonlinear_source_jacobian`).
"""
struct PhaseGrid{T<:Real}
    n_per_dim::Int
    phase_factors::Matrix{Complex{T}}
    jacobian_cosines::Matrix{T}   # cos(k_slot · ξ_point) / n_points
    jacobian_sines::Matrix{T}     # sin(k_slot · ξ_point) / n_points
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
    PhaseGrid(n_per_dim, phase_factors, real.(averaging_table), imag.(averaging_table),
              sum_slot, difference_slot, sum_conjugated, difference_conjugated)
end

n_points(grid::PhaseGrid) = size(grid.phase_factors, 1)

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
    component_profiles(quadrature, profiles)

Sum each component's channels over their polar harmonics, returning
`F[node, component] = F_q(x_node) = Σ_l R_{l,q} Ỹ_{l,|m_q|}(x_node)`.
"""
function component_profiles(quadrature::FieldQuadrature, profiles::AbstractVector)
    mode_set, harmonics = quadrature.mode_set, quadrature.angular.harmonics
    summed = zeros(promote_type(eltype(profiles), eltype(harmonics)), size(harmonics, 1), n_components(mode_set))
    for (component_index, channels) in enumerate(mode_set.component_channels)
        summed[:, component_index] = harmonics[:, channels] * profiles[channels]
    end
    return summed
end

"""
    reconstruct_field(quadrature, profiles)

The real field `Φ(x_node, ξ_point) = Σ_q F_q(x_node) e^{-i q·ξ_point} + c.c.` as a
matrix indexed `[point, node]` (phase point, polar node).
"""
reconstruct_field(quadrature::FieldQuadrature, profiles::AbstractVector) =
    2 .* real.(quadrature.phase.phase_factors * transpose(component_profiles(quadrature, profiles)))

"""
    project_phases(quadrature, grid_function)

Phase average `⟨e^{+i q·ξ} f⟩_ξ` of a grid function `f[point, node]`, returning
a matrix indexed `[component, node]`.
"""
project_phases(quadrature::FieldQuadrature, grid_function::AbstractMatrix) =
    (quadrature.phase.phase_factors' * grid_function) ./ n_points(quadrature.phase)

"""
    nonlinear_source(quadrature, potential, background, α, s, profiles)

The projected nonlinear source for every channel at one radial node,
`N_{l,q} = α² ∫ dx Ỹ_{l,|m_q|}(x) ρ²(r, x) ⟨e^{+i q·ξ} h(Φ)⟩_ξ`,
evaluated on the full reconstructed field without any expansion in the
self-interaction strength.
"""
function nonlinear_source(quadrature::FieldQuadrature, potential::AxionPotential,
                          background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    field = reconstruct_field(quadrature, profiles)
    projected = project_phases(quadrature, nonlinearity.(Ref(potential), field))
    return _project_polar(quadrature, background, s, projected) .* mass_coupling^2
end

nonlinear_source(quadrature::FieldQuadrature, ::FreePotential, ::KerrBackground, mass_coupling, s, profiles::AbstractVector) =
    zeros(promote_type(eltype(profiles), typeof(mass_coupling)), length(profiles))

"""Weighted polar projection `Σ_node w Ỹ_channel(x) ρ²(x) projected[component(channel), node]`."""
function _project_polar(quadrature::FieldQuadrature, background::KerrBackground, s, projected::AbstractMatrix)
    mode_set, grid = quadrature.mode_set, quadrature.angular
    weights = grid.weights .* rho2.(Ref(background), s, grid.x)
    return [sum(weights[node] * grid.harmonics[node, channel] * projected[mode_set.channel_component[channel], node]
                for node in eachindex(weights))
            for channel in 1:n_channels(mode_set)]
end

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

The coefficients at all polar nodes come from a single product of the phase
grid's averaging table with `h′(Φ)[point, node]` (two real BLAS products).
"""
function nonlinear_source_jacobian(quadrature::FieldQuadrature, potential::AxionPotential,
                                   background::KerrBackground, mass_coupling, s, profiles::AbstractVector)
    mode_set, grid, phase = quadrature.mode_set, quadrature.angular, quadrature.phase
    n_chan = n_channels(mode_set)
    T = promote_type(real(eltype(profiles)), typeof(mass_coupling), eltype(grid.x))
    wrt_real = zeros(Complex{T}, n_chan, n_chan)
    wrt_imag = zeros(Complex{T}, n_chan, n_chan)

    slope = nonlinearity_derivative.(Ref(potential), reconstruct_field(quadrature, profiles))  # h′(Φ)[point, node]
    fourier_coefficients = complex.(phase.jacobian_cosines * slope, phase.jacobian_sines * slope)  # ĥ′(k)[slot, node]
    polar_weights = grid.weights .* rho2.(Ref(background), s, grid.x) .* mass_coupling^2
    owner = mode_set.channel_component
    for node in eachindex(polar_weights), source in 1:n_chan, target in 1:n_chan
        target_component, source_component = owner[target], owner[source]
        plus = fourier_coefficients[phase.sum_slot[target_component, source_component], node]
        minus = fourier_coefficients[phase.difference_slot[target_component, source_component], node]
        phase.sum_conjugated[target_component, source_component] && (plus = conj(plus))
        phase.difference_conjugated[target_component, source_component] && (minus = conj(minus))
        prefactor = polar_weights[node] * grid.harmonics[node, target] * grid.harmonics[node, source]
        wrt_real[target, source] += prefactor * (plus + minus)
        wrt_imag[target, source] += prefactor * (-im) * (plus - minus)
    end
    return NonlinearSourceJacobian(wrt_real, wrt_imag)
end

function nonlinear_source_jacobian(quadrature::FieldQuadrature, ::FreePotential, ::KerrBackground,
                                   mass_coupling, s, profiles::AbstractVector)
    n_chan = n_channels(quadrature.mode_set)
    T = promote_type(real(eltype(profiles)), typeof(mass_coupling))
    NonlinearSourceJacobian(zeros(Complex{T}, n_chan, n_chan), zeros(Complex{T}, n_chan, n_chan))
end

"""
    phase_averaged_potential(quadrature, potential, profiles)

`⟨U(Φ)⟩_ξ` at each polar node for one radial node, as a vector over nodes.
"""
function phase_averaged_potential(quadrature::FieldQuadrature, potential::AxionPotential, profiles::AbstractVector)
    potential_on_grid = potential_energy.(Ref(potential), reconstruct_field(quadrature, profiles))
    return vec(sum(potential_on_grid; dims=1)) ./ n_points(quadrature.phase)
end
