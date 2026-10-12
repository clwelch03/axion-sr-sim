using FastGaussQuadrature: gausslegendre

"""
    polar_harmonic(l::Int, m::Int, cos_theta::Real)

Compute the normalized polar harmonic at `cos_theta` using a recurrence.

Require `0 ≤ m ≤ l` and `-1 ≤ cos_theta ≤ 1`. Include the Condon–Shortley
phase and normalize the squared harmonic to integrate to one over
`d(cos_theta)`. Build the `l = m` harmonic first, then increase `l`.

See also [`build_polar_values`](@ref).
"""
function polar_harmonic(l::Int, m::Int, cos_theta::Real)
    0 <= m <= l || throw(ArgumentError("require 0 <= m <= l"))
    -1 <= cos_theta <= 1 || throw(DomainError(cos_theta, "require -1 <= cos_theta <= 1"))

    current_harmonic = 1 / sqrt(2.0)  # Y_00
    for order in 1:m
        current_harmonic *= -sqrt((2 * order + 1) / (2 * order)) * sqrt(1 - cos_theta^2)
    end

    previous_harmonic = 0.0
    for degree in m:(l - 1)
        lower = degree == m ? 0.0 : sqrt((degree^2 - m^2) / (4 * degree^2 - 1))
        upper = sqrt(((degree + 1)^2 - m^2) / (4 * (degree + 1)^2 - 1))
        previous_harmonic, current_harmonic = current_harmonic, (cos_theta * current_harmonic - lower * previous_harmonic) / upper
    end
    return current_harmonic
end

"""
    build_polar_values(l_max, cosine_theta_nodes)

Tabulate polar harmonics for `0 ≤ m ≤ l ≤ l_max`.

Return a dictionary keyed by `(l, m)`, with vectors of values in
`cosine_theta_nodes` order. Use the normalization of [`polar_harmonic`](@ref).
"""
function build_polar_values(l_max, cosine_theta_nodes)
    return Dict(
        (l, m) => polar_harmonic.(l, m, cosine_theta_nodes)
        for l in 0:l_max for m in 0:l
    )
end

"""
    build_azimuthal_values(m_max, azimuth_nodes)

Tabulate `exp(im * m * ϕ)` at `azimuth_nodes` for `0 ≤ m ≤ m_max`.

Return a dictionary keyed by `m`, with vectors in node order. No additional
normalization factor is included; [`angular_grid`](@ref) supplies weights
for the normalized azimuthal measure `dϕ / (2π)`.
"""
function build_azimuthal_values(m_max, azimuth_nodes)
    return Dict(
        m => exp.(im * m .* azimuth_nodes)
        for m in 0:m_max
    )
end

"""
    angular_grid(n_polar, n_azimuth)

Construct a Gauss–Legendre polar grid and a uniform periodic azimuthal grid.

Return a named tuple with `cosine_theta_nodes`, `polar_weights`, `azimuth_nodes`,
and the scalar `azimuth_weight`. Polar weights integrate over `[-1, 1]` and
sum to two; azimuthal nodes span `[0, 2π)` with weight `1 / n_azimuth` each.

# Examples
```jldoctest
julia> grid = angular_grid(4, 8);

julia> (length(grid.cosine_theta_nodes), length(grid.azimuth_nodes))
(4, 8)

julia> grid.azimuth_weight
0.125
```
"""
function angular_grid(n_polar, n_azimuth)
    cosine_theta_nodes, polar_weights = gausslegendre(n_polar)
    azimuth_nodes = 2π .* (0:n_azimuth-1) ./ n_azimuth
    azimuth_weight = 1 / n_azimuth

    return (; cosine_theta_nodes, polar_weights,
              azimuth_nodes, azimuth_weight)
end
