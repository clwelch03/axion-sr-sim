"""
Bookkeeping for the multiphase expansion.

The field is a sum over Fourier components labelled by integer vectors `q`
(one entry per independently oscillating bound mode), and each component
carries a small set of retained polar indices `l`. We call one `(l, q)` pair a
*channel*: it is the unit that owns a radial profile `R_{l,q}`.
"""

"""
    BoundMode(n, l, m)

Label of an independently oscillating bound mode by its non-interacting
hydrogenic quantum numbers.
"""
struct BoundMode
    n::Int
    l::Int
    m::Int
    function BoundMode(n, l, m)
        n ≥ 1 && 0 ≤ l < n && abs(m) ≤ l || throw(ArgumentError("invalid |nlm⟩ = |$n$l$m⟩"))
        new(n, l, m)
    end
end

"""Hydrogenic real frequency `α(1 − α²/2n²)`, the small-`α` starting guess."""
hydrogenic_frequency(mode::BoundMode, mass_coupling) = mass_coupling * (1 - mass_coupling^2 / (2mode.n^2))

"""
    FourierComponent(q, m, ls)

One phase component `e^{-i q·ξ}` of the field: integer vector `q`, azimuthal
number `m = q·m_modes`, and retained polar indices `ls` (all with the parity
dictated by the modes that source the component).
"""
struct FourierComponent
    q::Vector{Int}
    m::Int
    ls::Vector{Int}
end

"""Whether `q` is a unit vector, i.e. one of the fundamental bound-mode components."""
is_fundamental(component::FourierComponent) = sum(abs, component.q) == 1

"""Fourier order `|q|₁` of a component."""
fourier_order(component::FourierComponent) = sum(abs, component.q)

"""
    ModeSet(modes; max_order=3, n_angular=1, even_potential=true)

Build the truncated set of retained Fourier components for `K` bound `modes`.

Components are all integer vectors with `0 < |q|₁ ≤ max_order` (odd `|q|₁` only
when the potential is even), keeping one representative of each `±q` pair with
`sum(q) > 0` so that `ω_q > 0`. Each component retains the `n_angular` lowest
polar indices `l ≥ |m_q|` of the correct parity, spaced by two because both
the Kerr `x²` coupling and the nonlinear source preserve reflection parity.
"""
struct ModeSet
    modes::Vector{BoundMode}
    components::Vector{FourierComponent}
    channel_component::Vector{Int}              # component index of each channel
    channel_l::Vector{Int}                      # polar index of each channel
    component_channels::Vector{UnitRange{Int}}  # channels owned by each component
end

function ModeSet(modes::AbstractVector{BoundMode}; max_order::Int=3, n_angular::Int=1,
                 even_potential::Bool=true)
    n_bound = length(modes)
    components = FourierComponent[]
    for q in _integer_vectors(n_bound, max_order)
        order = sum(abs, q)
        0 < order ≤ max_order || continue
        even_potential && iseven(order) && continue
        _is_representative(q) || continue
        azimuthal = sum(q[i] * modes[i].m for i in 1:n_bound)
        # Reflection parity of the component is the product of its sources' parities (−1)^(l−m).
        parity = mod(sum(q[i] * (modes[i].l - modes[i].m) for i in 1:n_bound), 2)
        lowest_l = abs(azimuthal) + parity
        push!(components, FourierComponent(q, azimuthal, [lowest_l + 2k for k in 0:n_angular-1]))
    end
    sort!(components; by=component -> (fourier_order(component), .-component.q))   # fundamentals first, in mode order

    channel_component = Int[]
    channel_l = Int[]
    component_channels = UnitRange{Int}[]
    for (component_index, component) in enumerate(components)
        first_channel = length(channel_l) + 1
        append!(channel_component, fill(component_index, length(component.ls)))
        append!(channel_l, component.ls)
        push!(component_channels, first_channel:length(channel_l))
    end
    ModeSet(collect(modes), components, channel_component, channel_l, component_channels)
end

"""All integer vectors of length `n` with entries in `−max_order:max_order`."""
_integer_vectors(n, max_order) =
    vec([collect(entries) for entries in Iterators.product(ntuple(_ -> -max_order:max_order, n)...)])

"""Pick the member of a `±q` pair with positive sum (or, if the sum vanishes, positive first nonzero entry)."""
function _is_representative(q)
    total = sum(q)
    total != 0 && return total > 0
    return q[findfirst(!iszero, q)] > 0
end

n_modes(mode_set::ModeSet) = length(mode_set.modes)
n_components(mode_set::ModeSet) = length(mode_set.components)
n_channels(mode_set::ModeSet) = length(mode_set.channel_l)

"""The component owning `channel`."""
channel_component(mode_set::ModeSet, channel) = mode_set.components[mode_set.channel_component[channel]]

"""Azimuthal number `m_q` of the component owning `channel`."""
channel_m(mode_set::ModeSet, channel) = channel_component(mode_set, channel).m

"""
    fundamental_channels(mode_set)

Channel index of the `(l_i, e_i)` profile of each bound mode `i`. These are the
channels whose leading outer coefficient is prescribed by the amplitude `A_i`.
"""
function fundamental_channels(mode_set::ModeSet)
    map(1:n_modes(mode_set)) do mode_index
        component_index = findfirst(component -> is_fundamental(component) && component.q[mode_index] == 1,
                                    mode_set.components)
        channels = mode_set.component_channels[component_index]
        channels[findfirst(==(mode_set.modes[mode_index].l), mode_set.channel_l[channels])]
    end
end

"""Complex frequency `ω̃_q = q·ω̃` of each component, given per-mode frequencies."""
component_frequencies(mode_set::ModeSet, mode_frequencies) =
    [sum(component.q .* mode_frequencies) for component in mode_set.components]

Base.show(io::IO, mode::BoundMode) = print(io, "|", mode.n, mode.l, mode.m, "⟩")
function Base.show(io::IO, mode_set::ModeSet)
    print(io, "ModeSet(", join(string.(mode_set.modes), ", "), "; ",
          n_components(mode_set), " components, ", n_channels(mode_set), " channels)")
end
