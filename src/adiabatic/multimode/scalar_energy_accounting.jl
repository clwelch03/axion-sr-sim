"""
    primary_boundary_fluxes(primary_branch, primary_energies)

Return `(; horizon, outer)`, vectors of the primary cloud's scalar energy
fluxes at each energy in `primary_energies`.

Supply `primary_branch` as a dictionary of saved column vectors with keys
`"energy"`, `"horizon_flux"`, and `"outer_flux"`. Its energies must be positive,
strictly increasing, and contain at least four points. Requested energies must
lie within that range; [`evaluate_spline`](@ref) does not extrapolate.

For each boundary b, interpolate F₁ᵇ/E₁ in the primary energy E₁ using
[`cubic_spline`](@ref), then multiply by the requested E₁. This is the same
small-field scaling used by [`two_cloud_coefficients`](@ref).

Fluxes are in code energy units per unit t/M, before physical normalization
or conversion to the saved growth-scaled time. `horizon` is positive into
the BH (negative for superradiant extraction); `outer` is positive outward
at the finite outer boundary, not a flux evaluated at infinity. These are
the primary's own fluxes, excluding the secondary and driven m = 0, 3 fields.
"""
function primary_boundary_fluxes(primary_branch, primary_energies)
    energy_nodes = primary_branch["energy"]

    function interpolate_flux(column)
        spline = cubic_spline(
            energy_nodes, primary_branch[column] ./ energy_nodes
        )
        return [
            energy * evaluate_spline(spline, energy).value
            for energy in primary_energies
        ]
    end

    return (
        horizon = interpolate_flux("horizon_flux"),
        outer = interpolate_flux("outer_flux"),
    )
end

"""
    secondary_boundary_fluxes(primary_branch, secondary_energy_table,
                              primary_energies, secondary_energies)

Return `(; horizon, outer)`, vectors of the secondary cloud's own scalar
energy fluxes at paired trajectory energies E₁, E₂ from `primary_energies`
and `secondary_energies`.

Supply dictionaries of saved column vectors: `primary_branch` requires
`"amplitude"` and `"energy"`; `secondary_energy_table` requires `"amplitude"`,
`"energy"`, `"horizon_flux"`, and `"outer_flux"`. The tables must describe the
same primary amplitudes in the same order. The secondary table's energy Ê₂
is for its reference amplitude, not the trajectory energy E₂.

For each boundary b, interpolate F₂,tableᵇ/Ê₂ in the primary energy using
[`cubic_spline`](@ref), then multiply by the actual secondary energy:

    F₂ᵇ(E₁, E₂) = E₂ * spline(F₂,tableᵇ/Ê₂)(E₁).

Energy and flux scale together with the squared reference amplitude, so this
ratio removes that arbitrary normalization. Supply positive tabulated Ê₂ and
at least four strictly increasing primary energy nodes. Requested primary
energies must lie within the nodes; no extrapolation is performed. Throw if
the trajectory lengths differ or the tables' amplitudes do not match.

Fluxes use the same code units and boundary signs as
[`primary_boundary_fluxes`](@ref): positive into the BH or outward at the
finite outer boundary. No physical or growth-time conversion is applied.
The driven m = 0, 3 fields are separate contributions and are not included.
"""
function secondary_boundary_fluxes(
    primary_branch, secondary_energy_table, primary_energies, secondary_energies
)
    length(primary_energies) == length(secondary_energies) ||
        throw(DimensionMismatch("require both cloud energies at each time"))
    all(isapprox.(primary_branch["amplitude"], secondary_energy_table["amplitude"];
                  rtol=1e-12, atol=0)) ||
        throw(ArgumentError("primary and secondary tables must have matching amplitudes"))

    function interpolate_flux(column)
        spline = cubic_spline(
            primary_branch["energy"],
            secondary_energy_table[column] ./ secondary_energy_table["energy"]
        )
        return [
            secondary_energy * evaluate_spline(spline, primary_energy).value
            for (primary_energy, secondary_energy) in
                zip(primary_energies, secondary_energies)
        ]
    end

    return (
        horizon = interpolate_flux("horizon_flux"),
        outer = interpolate_flux("outer_flux"),
    )
end

"""
    dissipative_boundary_fluxes(primary_branch, dissipative_table,
                                primary_energies, secondary_energies)

Return `(; horizon, outer)`, vectors of scalar energy fluxes from the driven
m = 0 field at paired trajectory energies E₁, E₂.

Supply saved column dictionaries: `primary_branch` requires `"amplitude"`
and `"energy"`; `dissipative_table` requires `"amplitude"`, `"secondary_energy"`,
`"horizon_flux"`, and `"outer_flux"`. Rows must match in primary amplitude.
The tabulated secondary energy Ê₂ corresponds to the reference amplitude
used to calculate the driven field.

For each boundary b, interpolate the coefficient F₀,tableᵇ/(E₁² Ê₂)
in primary energy using [`cubic_spline`](@ref), then reconstruct

    F₀ᵇ(E₁, E₂) = E₁² E₂ * spline(F₀,tableᵇ/(E₁² Ê₂))(E₁).

Factoring out E₁² gives a smooth dilute limit; the interpolated coefficient
retains the nonlinear dependence on E₁. Require positive tabulated energies,
at least four strictly increasing primary energy nodes, and requested
primary energies within their range. Throw if the trajectory lengths differ
or the tables' amplitudes do not match.

Fluxes use code energy units per unit t/M. No physical or growth-time conversion.
"""
function dissipative_boundary_fluxes(
    primary_branch, dissipative_table, primary_energies, secondary_energies
)
    length(primary_energies) == length(secondary_energies) ||
        throw(DimensionMismatch("require both cloud energies at each time"))
    all(isapprox.(primary_branch["amplitude"], dissipative_table["amplitude"];
                  rtol=1e-12, atol=0)) ||
        throw(ArgumentError("primary and dissipative tables must have matching amplitudes"))

    energy_nodes = primary_branch["energy"]
    normalization = energy_nodes .^ 2 .* dissipative_table["secondary_energy"]

    function interpolate_flux(column)
        spline = cubic_spline(energy_nodes, dissipative_table[column] ./ normalization)
        return [
            primary_energy^2 * secondary_energy *
                evaluate_spline(spline, primary_energy).value
            for (primary_energy, secondary_energy) in
                zip(primary_energies, secondary_energies)
        ]
    end

    return (
        horizon = interpolate_flux("horizon_flux"),
        outer = interpolate_flux("outer_flux"),
    )
end

"""
    radiative_boundary_fluxes(primary_branch, radiative_table,
                              primary_energies, secondary_energies)

Return `(; horizon, infinity)`, vectors of scalar energy fluxes from the
driven m = 3 field at paired trajectory energies E₁, E₂.

Supply saved column dictionaries: `primary_branch` requires `"amplitude"`
and `"energy"`; `radiative_table` requires `"amplitude"`, `"secondary_energy"`,
`"horizon_flux"`, and `"infinity_flux"`. Rows must match in primary amplitude.
The tabulated secondary energy Ê₂ corresponds to the reference amplitude
used to calculate the driven field.

For each boundary b, interpolate F₃,tableᵇ/(E₁ Ê₂²) in primary energy using
[`cubic_spline`](@ref), then reconstruct

    F₃ᵇ(E₁, E₂) = E₁ E₂² * spline(F₃,tableᵇ/(E₁ Ê₂²))(E₁).

Factoring out E₁ gives a smooth dilute limit; the interpolated coefficient
retains the nonlinear dependence on E₁. Require positive tabulated energies,
at least four strictly increasing primary energy nodes, and requested
primary energies within their range. Throw if the trajectory lengths differ
or the tables' amplitudes do not match.
"""
function radiative_boundary_fluxes(
    primary_branch, radiative_table, primary_energies, secondary_energies
)
    length(primary_energies) == length(secondary_energies) ||
        throw(DimensionMismatch("require both cloud energies at each time"))
    all(isapprox.(primary_branch["amplitude"], radiative_table["amplitude"];
                  rtol=1e-12, atol=0)) ||
        throw(ArgumentError("primary and radiative tables must have matching amplitudes"))

    energy_nodes = primary_branch["energy"]
    normalization = energy_nodes .* radiative_table["secondary_energy"] .^ 2

    function interpolate_flux(column)
        spline = cubic_spline(energy_nodes, radiative_table[column] ./ normalization)
        return [
            primary_energy * secondary_energy^2 *
                evaluate_spline(spline, primary_energy).value
            for (primary_energy, secondary_energy) in
                zip(primary_energies, secondary_energies)
        ]
    end

    return (
        horizon = interpolate_flux("horizon_flux"),
        infinity = interpolate_flux("infinity_flux"),
    )
end

"""
    cumulative_trapezoid(times, values)

Return the cumulative integral of sampled `values` over `times`, starting
at zero at the first sample. Both inputs must be nonempty vectors of equal
length, with finite entries and strictly increasing times. Return a vector
of Float64 values.
"""
function cumulative_trapezoid(times, values)
    length(times) == length(values) ||
        throw(DimensionMismatch("require one value per time"))
    isempty(times) && throw(ArgumentError("require at least one sample"))
    all(isfinite, times) && all(isfinite, values) ||
        throw(ArgumentError("times and values must be finite"))
    all(diff(times) .> 0) ||
        throw(ArgumentError("times must increase strictly"))

    integral = zeros(Float64, length(times))
    for index in 2:length(times)
        spacing = times[index] - times[index - 1]
        integral[index] = integral[index - 1] +
            spacing * (values[index - 1] + values[index]) / 2
    end
    return integral
end

"""
    scalar_boundary_torques(coefficients, primary_energies,
                            primary_fluxes, secondary_fluxes, radiative_fluxes)

Return `(; horizon, outward)`, the total scalar angular-momentum fluxes
along a trajectory. Supply coefficient splines from
[`two_cloud_coefficients`](@ref) and boundary energy-flux vectors from the
primary, secondary, and radiative reconstruction functions, all in the same
trajectory order.

Use F_J = (m/omega) F_E with real code frequencies: 1/omega₁ for the
primary, 2/omega₂ for the secondary, and 3/(2omega₂ - omega₁) for the
driven m = 3 field. The driven m = 0 field carries zero torque.

The returned rates are in code angular-momentum units per unit t/M.
Horizon torque is positive into the BH. Outward torque includes finite
outer-boundary currents for the bound clouds and the m = 3 infinity flux.
"""
function scalar_boundary_torques(
    coefficients, primary_energies, primary_fluxes, secondary_fluxes, radiative_fluxes
)
    primary_omega = [
        evaluate_spline(coefficients.primary_omega, energy).value
        for energy in primary_energies
    ]
    secondary_omega = [
        evaluate_spline(coefficients.secondary_omega, energy).value
        for energy in primary_energies
    ]
    radiative_omega = 2secondary_omega .- primary_omega

    return (
        horizon = primary_fluxes.horizon ./ primary_omega .+
            2 .* secondary_fluxes.horizon ./ secondary_omega .+
            3 .* radiative_fluxes.horizon ./ radiative_omega,
        outward = primary_fluxes.outer ./ primary_omega .+
            2 .* secondary_fluxes.outer ./ secondary_omega .+
            3 .* radiative_fluxes.infinity ./ radiative_omega,
    )
end