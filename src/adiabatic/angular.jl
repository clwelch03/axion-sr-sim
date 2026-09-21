using FastGaussQuadrature: gausslegendre

"""
    polar_harmonic(l, m, cos_theta)

Compute the polar part of the spherical harmonic Y_{lm}(cos theta) via the recurrence relation.

Builds Y_{mm} first, then increases l until we reach the desired harmonic. Normalized over d(cos theta).

Return
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

Construct a grid of polar harmonic values for l <= l_max, m <= l.
"""
function build_polar_values(l_max, cosine_theta_nodes)
    return Dict(
        (l, m) => polar_harmonic.(l, m, cosine_theta_nodes)
        for l in 0:l_max for m in 0:l
    )
end

"""
    build_azimuthal_values(m_max, azimuth_nodes)

Construct a grid of azimuthal harmonic values for m <= m_max.
"""
function build_azimuthal_values(m_max, azimuth_nodes)
    return Dict(
        m => exp.(im * m .* azimuth_nodes)
        for m in 0:m_max
    )
end

function angular_grid(n_polar, n_azimuth)
    cosine_theta_nodes, polar_weights = gausslegendre(n_polar)
    azimuth_nodes = 2π .* (0:n_azimuth-1) ./ n_azimuth
    azimuth_weight = 1 / n_azimuth

    return (; cosine_theta_nodes, polar_weights,
              azimuth_nodes, azimuth_weight)
end