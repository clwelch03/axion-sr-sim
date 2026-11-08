"""
Polar spherical harmonics `Ỹ_{lm}(x)`, `x = cos θ`, normalised so that
`∫₋₁¹ Ỹ_{lm}² dx = 1`. The azimuthal factor `e^{imφ}` is carried separately by
each Fourier component.
"""

"""
    polar_harmonic(l, m, x)

Normalised associated Legendre function `Ỹ_{lm}(x)` for `0 ≤ m ≤ l`, computed
with the standard stable upward recurrence in `l` (Condon–Shortley phase included).
"""
function polar_harmonic(l::Int, m::Int, x)
    0 ≤ m ≤ l || throw(ArgumentError("need 0 ≤ m ≤ l, got l=$l m=$m"))
    # Seed with the normalised diagonal term Ỹ_mm, then climb in l.
    diagonal = sqrt(one(x) / 2)
    for k in 1:m
        diagonal *= -sqrt((2k + 1) / (2k)) * sqrt(1 - x^2)
    end
    l == m && return diagonal
    previous, current = diagonal, x * sqrt(2m + 3) * diagonal
    for k in (m + 2):l
        forward = sqrt((4k^2 - 1) / (k^2 - m^2))
        backward = sqrt(((k - 1)^2 - m^2) / (4(k - 1)^2 - 1))
        previous, current = current, forward * (x * current - backward * previous)
    end
    return current
end

"""
    polar_harmonic_dtheta(l, m, x)

`∂Ỹ_{lm}/∂θ` expressed through `x = cos θ`, using the identity
`(1 − x²) dỸ/dx = −l x Ỹ_{lm} + √((2l+1)/(2l−1)) √((l−m)(l+m)) Ỹ_{l−1,m}`
and `∂/∂θ = −sin θ ∂/∂x`.
"""
function polar_harmonic_dtheta(l::Int, m::Int, x)
    sin_theta = sqrt(1 - x^2)
    harmonic = polar_harmonic(l, m, x)
    lower_harmonic = l > m ? sqrt((2l + 1) / (2l - 1) * (l - m) * (l + m)) * polar_harmonic(l - 1, m, x) : zero(harmonic)
    sin2_times_dharmonic_dx = -l * x * harmonic + lower_harmonic
    return -sin2_times_dharmonic_dx / sin_theta
end

"""
    AngularGrid(n_nodes, mode_set)

Gauss–Legendre quadrature in `x = cos θ` together with tables of `Ỹ_{l,|m_q|}`
and `∂_θ Ỹ_{l,|m_q|}` evaluated for every channel. `harmonics[node, channel]`.
"""
struct AngularGrid{T<:Real}
    x::Vector{T}
    weights::Vector{T}
    harmonics::Matrix{T}
    harmonics_dtheta::Matrix{T}
end

function AngularGrid(n_nodes::Int, mode_set::ModeSet)
    x, weights = gausslegendre(n_nodes)
    channels = 1:n_channels(mode_set)
    harmonics = [polar_harmonic(mode_set.channel_l[channel], abs(channel_m(mode_set, channel)), node)
                 for node in x, channel in channels]
    harmonics_dtheta = [polar_harmonic_dtheta(mode_set.channel_l[channel], abs(channel_m(mode_set, channel)), node)
                        for node in x, channel in channels]
    AngularGrid(x, weights, harmonics, harmonics_dtheta)
end

"""
    mixing_matrix(grid, mode_set, component_index)

The matrix `C^m_{ll'} = ∫ Ỹ_{l,|m|} x² Ỹ_{l',|m|} dx` over the channels of one
component, which couples its polar indices through the Kerr `a² x²` term.
"""
function mixing_matrix(grid::AngularGrid, mode_set::ModeSet, component_index::Int)
    channels = mode_set.component_channels[component_index]
    harmonics = @view grid.harmonics[:, channels]
    return harmonics' * (grid.weights .* grid.x .^ 2 .* harmonics)
end
