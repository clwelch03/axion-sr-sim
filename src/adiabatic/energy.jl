"""
    energy_per_tortoise_radius(state, tortoise_radius, frequency, spin, alpha,
                               modes, polar_values, azimuthal_values, grid;
                               fundamental_m=1)

Compute the angular-integrated energy per unit tortoise radius at time zero.


Parallels Eq. (29) from OTT. There are five terms in (29), which are evaluated
in this function as follows:

* (D_t phi)^2       -> `sampled_energy`
* (D_varphi phi)^2  -> `angular_energy` and `rotation_energy`
* (D_r phi)^2       -> `radial_energy`
* (D_theta phi)^2   -> `angular_energy`
* V(phi)            -> `sampled_energy`

Supply `state = [R; Q]`, with unrescaled radial profiles `R` and slopes
`Q = dR/dr_*`, both in `modes` order. Each mode is a unique `(l, m)` pair
with positive `m`. Use units `M = 1`, with `spin = a/M`, `alpha = M*μ`,
and complex fundamental `frequency`; each harmonic has frequency
`(m/fundamental_m)*frequency`.

Supply angular tables from [`build_polar_values`](@ref) and
[`build_azimuthal_values`](@ref) on `grid` from [`angular_grid`](@ref).
Use the cosine potential `V(φ) = alpha^2*(1-cos(φ))`.

Evaluate time and potential terms by angular quadrature; integrate the
remaining derivative terms using harmonic orthogonality. This implements
the explicit energy terms in PDF Eq. (29), with the measure from Eq. (25).
Radial integration is left to the caller.

Return the real scalar `d(E/Fa²)/dr_*`, including the radial Jacobian
`dr/dr_* = Δ/(r²+a²)` and the full azimuthal measure `dϕ`.

# Arguments
- `fundamental_m=1`: positive integer azimuthal number of the fundamental.

See also [`reconstruct_field`](@ref), [`tortoise_geometry`](@ref).
"""
function energy_per_tortoise_radius(
    state, tortoise_radius, frequency, spin, alpha, modes,
    polar_values, azimuthal_values, grid;
    fundamental_m=1,
)
    mode_count = length(modes)
    length(state) == 2mode_count ||
        throw(DimensionMismatch("require state = [R; Q]"))

    geometry = tortoise_geometry(tortoise_radius, spin)
    radial_factor = geometry.radial_factor       # B = r²+a²
    jacobian = geometry.dr_drstar                # Δ/B

    radial_values = Dict(
        mode => state[index] for (index, mode) in enumerate(modes)
    )

    # Differentiate exp(-i*n*omega*t), where n=m/fundamental_m.
    time_values = Dict(
        (l, m) => -im * (m / fundamental_m) * frequency *
                  radial_values[(l, m)]
        for (l, m) in modes
    )


    field = reconstruct_field(
        radial_values, polar_values, azimuthal_values, grid,
    )
    time_derivative = reconstruct_field(
        time_values, polar_values, azimuthal_values, grid,
    )

    sine_squared = 1 .- grid.cosine_theta_nodes.^2
    rho_squared =
        geometry.radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2

    # Time-derivative and potential terms require angular quadrature.
    time_coefficient =
        (radial_factor .- jacobian .* spin^2 .* sine_squared) ./ 2

    # Equivalent to alpha²*(1-cos(field)), without small-field cancellation.
    potential = 2 .* alpha^2 .* sin.(field ./ 2).^2

    sampled_density =
        time_coefficient .* time_derivative.^2 .+
        jacobian .* rho_squared .* potential

    # Existing azimuthal weights integrate dvarphi/(2pi).
    sampled_energy = 2π * grid.azimuth_weight *
        sum(grid.polar_weights .* sampled_density)

    # Remaining angular integrals follow from harmonic orthogonality.
    radial_energy = 2π * radial_factor *
        sum(abs2, state[mode_count+1:end])

    angular_energy = 2π * jacobian * sum(
        l * (l + 1) * abs2(radial_values[(l, m)])
        for (l, m) in modes
    )

    rotation_energy = -2π * spin^2 / radial_factor * sum(
        m^2 * abs2(radial_values[(l, m)])
        for (l, m) in modes
    )

    return sampled_energy + radial_energy + angular_energy + rotation_energy
end

"""
    integrate_energy_branch(branch, frequency, spin, alpha, modes,
                            polar_values, azimuthal_values, grid;
                            fundamental_m=1)

Integrate the energy density over one retained radial branch.

Supply `branch.positions` and `branch.states` from
[`integrate_radial`](@ref) with `save_profile=true`. Each state column
contains `[R; Q]` in `modes` order. Accept strictly increasing or decreasing
positions; orient the integral toward increasing tortoise radius.

Return the real scalar `E/Fa²` over the branch's finite radial interval,
using composite Simpson quadrature. Require uniform spacing and an even
number of intervals on each branch. This quadrature is fourth order for
smooth energy densities; check radial spacing and endpoint dependence
separately.

# Arguments
- `fundamental_m=1`: positive integer azimuthal number of the fundamental.

See also [`energy_per_tortoise_radius`](@ref).
"""
function integrate_energy_branch(
    branch, frequency, spin, alpha, modes,
    polar_values, azimuthal_values, grid;
    fundamental_m=1,
)
    positions = branch.positions
    length(positions) >= 2 ||
        throw(ArgumentError("require at least two radial samples"))
    size(branch.states, 2) == length(positions) ||
        throw(DimensionMismatch("require one state column per position"))

    spacings = diff(positions)
    if all(spacings .> 0)
        order = eachindex(positions)
    elseif all(spacings .< 0)
        order = reverse(eachindex(positions))
    else
        throw(ArgumentError("positions must be strictly monotonic"))
    end

    ordered_positions = positions[order]
    densities = [
        energy_per_tortoise_radius(
            view(branch.states, :, index), positions[index],
            frequency, spin, alpha, modes,
            polar_values, azimuthal_values, grid;
            fundamental_m,
        )
        for index in order
    ]

    # Simpson's rule
    interval_count = length(ordered_positions) - 1
    iseven(interval_count) ||
        throw(ArgumentError("Simpson's rule requires an even number of intervals"))

    spacing = (last(ordered_positions) - first(ordered_positions)) / interval_count
    all(
        step -> isapprox(step, spacing; rtol=1e-10, atol=0),
        diff(ordered_positions),
    ) || throw(ArgumentError("Simpson's rule requires uniformly spaced positions"))

    # Julia index 2 corresponds to the first interior point, with weight 4.
    weighted_density =
        first(densities) + last(densities) +
        4 * sum(densities[2:2:end-1]) +
        2 * sum(densities[3:2:end-2])

    return (spacing / 3) * weighted_density
end
