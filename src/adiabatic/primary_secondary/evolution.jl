"""
    cubic_spline(nodes, values)

Return the not-a-knot cubic spline through the points `(nodes[i], values[i])`
as a named tuple `(; nodes, values, second_derivatives)`, for
[`evaluate_spline`](@ref).

On each interval the spline is the cubic with the given end values and second
derivatives `M_i`. The `M_i` follow from a linear system: at interior nodes the
first derivative is continuous; at the second and second-to-last nodes the
third derivative is also continuous ("not-a-knot"), which replaces the free
end conditions. With that choice the spline reproduces any cubic exactly.
Nodes must increase strictly; at least four are required.
"""
function cubic_spline(nodes, values)
    count = length(nodes)
    count == length(values) || throw(DimensionMismatch("require one value per node"))
    count >= 4 || throw(ArgumentError("require at least four nodes"))
    spacings = diff(nodes)
    all(spacings .> 0) || throw(ArgumentError("nodes must increase strictly"))
    slopes = diff(values) ./ spacings

    matrix = zeros(count, count)
    right_side = zeros(count)
    # Continuity of the first derivative at interior node i.
    for i in 2:count-1
        matrix[i, i-1] = spacings[i-1]
        matrix[i, i] = 2(spacings[i-1] + spacings[i])
        matrix[i, i+1] = spacings[i]
        right_side[i] = 6(slopes[i] - slopes[i-1])
    end
    # Not-a-knot: (M2 - M1)/h1 = (M3 - M2)/h2, and the same at the other end.
    matrix[1, 1:3] = [spacings[2], -(spacings[1] + spacings[2]), spacings[1]]
    matrix[count, count-2:count] =
        [spacings[end], -(spacings[end-1] + spacings[end]), spacings[end-1]]

    return (; nodes, values, second_derivatives=matrix \ right_side)
end

"""
    evaluate_spline(spline, position)

Return `(; value, derivative)` of a [`cubic_spline`](@ref) at `position`, which
must lie within the nodes (no extrapolation).
"""
function evaluate_spline(spline, position)
    (; nodes, values, second_derivatives) = spline
    first(nodes) <= position <= last(nodes) ||
        throw(DomainError(position, "outside the spline's nodes"))
    i = clamp(searchsortedlast(nodes, position), 1, length(nodes) - 1)
    spacing = nodes[i+1] - nodes[i]

    # Weights of the left and right node; each runs from 1 to 0 across the interval.
    left_weight = (nodes[i+1] - position) / spacing
    right_weight = (position - nodes[i]) / spacing
    value = left_weight * values[i] + right_weight * values[i+1] +
            ((left_weight^3 - left_weight) * second_derivatives[i] +
             (right_weight^3 - right_weight) * second_derivatives[i+1]) * spacing^2 / 6
    derivative = (values[i+1] - values[i]) / spacing +
                 (-(3left_weight^2 - 1) * second_derivatives[i] +
                  (3right_weight^2 - 1) * second_derivatives[i+1]) * spacing / 6
    return (; value, derivative)
end

const TWO_CLOUD_COEFFICIENT_NAMES = (
    :primary_omega, :secondary_omega, :primary_flux_per_energy,
    :secondary_flux_per_energy, :m0_flux_coefficient, :m3_flux_coefficient,
)

"""
    two_cloud_coefficients(table; points=eachindex(table.primary_energy))

Return a named tuple of [`cubic_spline`](@ref)s in the primary energy E₁, one
per coefficient of the balance equations, built from the saved `points` of
`table` (the vectors returned by `read_two_cloud_table`):

- `primary_omega`, `secondary_omega`: Re ω₁, Re ω₂;
- `primary_flux_per_energy`: the primary's own flux F₁/E₁;
- `secondary_flux_per_energy`: the secondary's own flux per unit E₂;
- `m0_flux_coefficient`: the m = 0 flux per unit E₂ E₁²;
- `m3_flux_coefficient`: the m = 3 flux per unit E₂² E₁.

Each flux is divided by its small-field power, so it tends to a constant as
E₁ → 0; the coefficients are smooth in E₁ ∝ A₁², which is why the splines
use E₁ and not ln E₁. Pass `points` to build from a subset, for checks.
"""
function two_cloud_coefficients(table; points=eachindex(table.primary_energy))
    nodes = table.primary_energy[points]
    return NamedTuple{TWO_CLOUD_COEFFICIENT_NAMES}(Tuple(
        cubic_spline(nodes, getproperty(table, name)[points])
        for name in TWO_CLOUD_COEFFICIENT_NAMES
    ))
end

"""
    two_cloud_rates(coefficients, primary_energy, secondary_energy;
                    channels=(:primary, :secondary, :m0, :m3))

Return `(; primary_rate, secondary_rate, energy_flux, angular_momentum_flux)`:
dE₁/dt and dE₂/dt from energy and angular-momentum balance (Omiya et al.
2211.01949 Eqs. (62)-(63), written in the cloud energies; see
docs/multimode-adiabatic-design.md, "Piece 6 design"), and the total fluxes.

With J = J₁(E₁) + 2E₂/ω₂(E₁) and dJ₁ = dE₁/ω₁ (Eq. (72)), the balances are

    dE₁/dt + dE₂/dt = -F_E,
    (1/ω₁ - (2E₂/ω₂²) dω₂/dE₁) dE₁/dt + (2/ω₂) dE₂/dt = -F_J.

F_E sums the energy fluxes of the selected `channels`: the primary's own
(∝ E₁), the secondary's own (∝ E₂), the m = 0 field (∝ E₁² E₂), and the m = 3
field (∝ E₁ E₂²). In F_J each carries m/ω: 1/ω₁, 2/ω₂, 0, 3/ω₃, with
ω₃ = 2ω₂ - ω₁. Positive flux means energy lost. Times are in units of M.
"""
function two_cloud_rates(coefficients, primary_energy, secondary_energy;
                         channels=(:primary, :secondary, :m0, :m3))
    value(name) = evaluate_spline(getproperty(coefficients, name), primary_energy).value
    primary_omega = value(:primary_omega)
    secondary = evaluate_spline(coefficients.secondary_omega, primary_energy)
    secondary_omega = secondary.value
    radiative_omega = 2secondary_omega - primary_omega

    # Energy flux of each channel, and its m/ω.
    fluxes = (
        primary=(value(:primary_flux_per_energy) * primary_energy, 1 / primary_omega),
        secondary=(value(:secondary_flux_per_energy) * secondary_energy, 2 / secondary_omega),
        m0=(value(:m0_flux_coefficient) * primary_energy^2 * secondary_energy, 0.0),
        m3=(value(:m3_flux_coefficient) * primary_energy * secondary_energy^2,
            3 / radiative_omega),
    )
    all(channel -> haskey(fluxes, channel), channels) ||
        throw(ArgumentError("channels must be among $(keys(fluxes))"))
    energy_flux = sum(first(fluxes[channel]) for channel in channels)
    angular_momentum_flux = sum(prod(fluxes[channel]) for channel in channels)

    # dJ/dE₁ at fixed E₂ includes the shift of J₂ = 2E₂/ω₂ with ω₂(E₁).
    balance = [1.0 1.0;
               1 / primary_omega - 2secondary_energy * secondary.derivative / secondary_omega^2 2 / secondary_omega]
    primary_rate, secondary_rate = balance \ [-energy_flux, -angular_momentum_flux]
    return (; primary_rate, secondary_rate, energy_flux, angular_momentum_flux)
end

"""
    evolve_two_clouds(coefficients, primary_energy, secondary_energy;
                      time_unit, final_time, time_step, channels=...)

Integrate [`two_cloud_rates`](@ref) from the initial energies with RK4 in
(ln E₁, ln E₂), which follow the clouds' exponential growth and decay smoothly.
Time is measured in units of `time_unit` (a rate in 1/M, e.g. the linear
growth rate Im ω₁), up to `final_time` with step `time_step`.

Return `(; times, primary_energies, secondary_energies)`. Throw if E₁ leaves
the range of the saved branch.
"""
function evolve_two_clouds(coefficients, primary_energy, secondary_energy;
                           time_unit, final_time, time_step,
                           channels=(:primary, :secondary, :m0, :m3))
    count = max(1, round(Int, final_time / time_step))
    spacing = final_time / count
    # d(ln E)/dτ with τ = time_unit * t.
    function log_rates(state)
        energies = exp.(state)
        rates = two_cloud_rates(coefficients, energies...; channels)
        return [rates.primary_rate, rates.secondary_rate] ./ energies ./ time_unit
    end

    state = log.([primary_energy, secondary_energy])
    times = collect(range(0, final_time; length=count + 1))
    states = Matrix{Float64}(undef, 2, count + 1)
    states[:, 1] = state
    for index in 1:count
        k1 = log_rates(state)
        k2 = log_rates(state .+ (spacing / 2) .* k1)
        k3 = log_rates(state .+ (spacing / 2) .* k2)
        k4 = log_rates(state .+ spacing .* k3)
        state = state .+ (spacing / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        states[:, index+1] = state
    end
    return (; times, primary_energies=exp.(states[1, :]), secondary_energies=exp.(states[2, :]))
end
