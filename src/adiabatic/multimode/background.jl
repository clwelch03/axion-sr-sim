"""
    azimuthal_average(field, grid)

Average an angular field over the azimuthal nodes of `grid`.

Supply `field` as a matrix with polar rows and azimuthal columns, as returned
by [`reconstruct_field`](@ref). Return a vector with one entry per polar node.
The uniform azimuthal rule integrates `exp(im*m*φ)` exactly for
`0 ≤ abs(m) < n_azimuth`.
"""
azimuthal_average(field, grid) = vec(sum(field; dims=2)) .* grid.azimuth_weight

"""
    averaged_mass_deficit(radial_values, polar_values, azimuthal_values, grid)

Return `<1 - cos ϕ₁>` at each polar node, at one radius: the average over
one period of the primary, computed as the azimuthal average at t = 0.

`radial_values` holds the primary's `R_lm` at that radius, keyed by `(l, m)`,
in the convention of [`reconstruct_field`](@ref). The secondary of
Omiya et al. 2211.01949 Eq. (61) sees the effective mass squared
`mu^2 * (1 - <1 - cos ϕ₁>)`.
"""
function averaged_mass_deficit(radial_values, polar_values, azimuthal_values, grid)
    field = reconstruct_field(radial_values, polar_values, azimuthal_values, grid)

    # 2 sin²(ϕ/2) equals 1 - cos ϕ without cancellation when ϕ is small.
    mass_deficit = 2 .* sin.(field ./ 2) .^ 2

    return azimuthal_average(mass_deficit, grid)
end

"""
    background_coupling_matrix(mass_deficit_average, radius, spin, modes,
                               polar_values, grid)

Return `V` with `V[i, j] = ∫ dx Y_i(x) Y_j(x) (r^2 + a^2 x^2) c0(x)` at one radius,
for the secondary's `modes` (all with one azimuthal number).

Supply `mass_deficit_average` as `c0` at the polar nodes of `grid`, from
[`averaged_mass_deficit`](@ref). `c0` is axisymmetric, so it cannot couple
different `m`. In the secondary's radial equation the background enters as
`alpha^2 * V * R` in the slot of the nonlinear source `N`
(Omiya et al. 2211.01949 Eq. (61), including the factor mu^2).
"""
function background_coupling_matrix(mass_deficit_average, radius, spin, modes,
                                    polar_values, grid)
    azimuthal_number = last(first(modes))
    all(mode -> last(mode) == azimuthal_number, modes) ||
        throw(ArgumentError("all modes must share one azimuthal number"))

    # Polar weight times rho^2 = r^2 + a^2 x^2 times c0, at each polar node.
    weighted_deficit = grid.polar_weights .*
        (radius^2 .+ spin^2 .* grid.cosine_theta_nodes .^ 2) .* mass_deficit_average

    # Rows follow the output mode, columns the input mode.
    return [
        sum(weighted_deficit .* polar_values[mode] .* polar_values[other_mode])
        for mode in modes, other_mode in modes
    ]
end

"""
    mass_deficit_harmonic(radial_values, polar_values, azimuthal_values, grid,
                          azimuthal_number)

Return the coefficient of `exp(im * azimuthal_number * φ)` in `1 - cos ϕ₁` at
each polar node, at one radius and time zero.

For the primary (m₁ = 1), `exp(+im*k*φ)` at t = 0 is `exp(-im*k*ψ₁)`, with
ψ₁ = ω₁t - φ, so `azimuthal_number = k` gives the part of `1 - cos ϕ₁`
varying as `exp(-im*k*ψ₁)`. `azimuthal_number = 0` is
[`averaged_mass_deficit`](@ref). The m = 0 field at 2ω₁ - ω₂ is driven by
`azimuthal_number = 2` (Omiya et al. 2211.01949 Eqs. (51)-(52)).
"""
function mass_deficit_harmonic(radial_values, polar_values, azimuthal_values, grid,
                               azimuthal_number)
    field = reconstruct_field(radial_values, polar_values, azimuthal_values, grid)

    # 2 sin²(ϕ/2) equals 1 - cos ϕ without cancellation when ϕ is small.
    mass_deficit = 2 .* sin.(field ./ 2) .^ 2

    # Averaging against conj(exp(im*m*φ)) keeps the exp(+im*m*φ) term.
    return azimuthal_average(
        mass_deficit .* transpose(conj.(azimuthal_values[azimuthal_number])), grid,
    )
end
