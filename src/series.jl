"""
Linear asymptotic expansions of the radial profiles at the two ends of the
domain. Both are built the same way: the radial equation multiplied through by
`Δ` is applied to a truncated series using small polynomial arithmetic, and the
coefficients are determined order by order with `_solve_series`. Writing the
operators rather than hand-derived recursions keeps the expansions transparent
and lets the same code verify itself (the leading residual must vanish by the
choice of the leading exponent).
"""

# --- tiny dense polynomial helpers (ascending coefficients) ------------------

_polymul(left::AbstractVector, right::AbstractVector) =
    [sum(left[i] * right[k - i + 1] for i in max(1, k - length(right) + 1):min(k, length(left)))
     for k in 1:(length(left) + length(right) - 1)]

"""Sum of two polynomials of possibly different lengths (ascending coefficients)."""
function _polyadd(left::AbstractVector, right::AbstractVector)
    result = zeros(promote_type(eltype(left), eltype(right)), max(length(left), length(right)))
    result[1:length(left)] .+= left
    result[1:length(right)] .+= right
    return result
end

# --- shared order-by-order solver -------------------------------------------

"""
    SeriesTerm(polynomial, shift, order_factors, channel_map)

One term of a linear operator on a series `Σ_j c_j x^j` whose coefficients `c_j` are vectors over channels:

    polynomial(x) · x^shift · Σ_j order_factors[j + 1] (channel_map c_j) x^j.

The operators of the horizon and outer expansions are sums of a few such terms (`_horizon_terms`,
`_outer_terms`), so one coefficient of the result costs `O(length(polynomial) n_channels²)` per term.
"""
struct SeriesTerm{T<:Number}
    polynomial::Vector{T}
    shift::Int
    order_factors::Vector{T}
    channel_map::Matrix{T}
end

"""Coefficient of `x^row_order` of `term` applied to the series rows `0:known_order` of `coefficients`, added to `result`."""
function _add_term_row!(result::AbstractVector, term::SeriesTerm, coefficients::AbstractMatrix, row_order::Int, known_order::Int)
    for (index, polynomial_coefficient) in enumerate(term.polynomial)
        order = row_order - term.shift - (index - 1)
        0 ≤ order ≤ known_order || continue
        factor = polynomial_coefficient * term.order_factors[order + 1]
        for target in eachindex(result), source in eachindex(result)
            result[target] += factor * term.channel_map[target, source] * coefficients[order + 1, source]
        end
    end
    return result
end

"""Matrix of the coefficient of `x^row_order` of `term` with respect to the series row `c_order`, added to `coupling`."""
function _add_term_coupling!(coupling::AbstractMatrix, term::SeriesTerm, row_order::Int, order::Int)
    index = row_order - term.shift - order + 1
    1 ≤ index ≤ length(term.polynomial) || return coupling
    coupling .+= (term.polynomial[index] * term.order_factors[order + 1]) .* term.channel_map
    return coupling
end

"""
    _solve_series(terms, leading, order, row_shift)

Fill a coefficient matrix `coefficients[j+1, channel]`, `j = 0:order`, such
that the operator `Σ terms` (`SeriesTerm`) applied to the series vanishes
through order `order + row_shift`. Unknown row `j` is determined by the
operator's coefficient `j + row_shift`, which is linear in that row with the
small channel-coupling matrix from `_add_term_coupling!`.
"""
function _solve_series(terms::AbstractVector{<:SeriesTerm}, leading::AbstractVector, order::Int, row_shift::Int)
    n_chan = length(leading)
    T = promote_type(eltype(leading), eltype(first(terms).polynomial))
    coefficients = zeros(T, order + 1, n_chan)
    coefficients[1, :] = leading
    residual, coupling = zeros(T, n_chan), zeros(T, n_chan, n_chan)
    for j in 1:order
        fill!(residual, 0)
        fill!(coupling, 0)
        for term in terms
            _add_term_row!(residual, term, coefficients, j + row_shift, j - 1)
            _add_term_coupling!(coupling, term, j + row_shift, j)
        end
        coefficients[j + 1, :] = -(coupling \ residual)
    end
    return coefficients
end

"""
Channel map `diag(constants) + coupling_strength · mixing` of the terms of the radial equation that couple channels:
`ΔB` has a channel-dependent constant and the angular coupling is `a²(ω²−α²) C`.
"""
function _channel_map(T::Type, constants::AbstractVector, coupling_strength, mixing::AbstractMatrix)
    channel_map = Matrix{T}(coupling_strength .* mixing)
    for channel in eachindex(constants)
        channel_map[channel, channel] += constants[channel]
    end
    return channel_map
end

# --- horizon ------------------------------------------------------------------

"""
    HorizonSeries

Frobenius expansion `R_{l,q}(s) = s^p Σ_j h^{(j)}_l s^j` of one component's
channels near the horizon, with `p = -i (r₊²+a²)(ω_q − m_q Ω_H)/(r₊ − r₋)`
selecting the ingoing branch. `coefficients[j+1, channel]`.
"""
struct HorizonSeries{P<:Number, T<:Number}
    horizon_power::P
    coefficients::Matrix{T}
end

"""Ingoing-branch horizon exponent `p_q` for frequency `ω` and azimuthal number `m`."""
horizon_exponent(background::KerrBackground, frequency, m) =
    -im * (2background.rplus / background.horizon_gap) * (frequency - m * background.horizon_angular_velocity)

"""
    HorizonSeries(background, α, ω, m, ls, mixing, leading, order)

Build the horizon expansion for a component with frequency `ω`, azimuthal
number `m`, polar indices `ls`, mixing matrix `C`, and leading coefficients.
"""
function HorizonSeries(background::KerrBackground, mass_coupling, frequency, m::Int, ls::AbstractVector{Int},
                       mixing::AbstractMatrix, leading::AbstractVector, order::Int)
    p = horizon_exponent(background, frequency, m)
    terms = _horizon_terms(background, mass_coupling, frequency, m, ls, mixing, p, order)
    HorizonSeries(p, _solve_series(terms, leading, order, 0))
end

"""
`SeriesTerm`s of `Δ(ΔR')' + ΔB R + Δ a²(ω²−α²) C R` applied to `R = s^p Σ_j c_j s^j`, as coefficients of
`s^{p+J}`. With `Δ = s(s + gap)` every coefficient is a polynomial in `s`:

    Δ²R'' + ΔΔ'R' = (gap² + 2gap s + s²) Σ (p+j)(p+j−1) c_j s^{p+j} + (gap² + 3gap s + 2s²) Σ (p+j) c_j s^{p+j},
    ΔB = g² + Δ (2amω − a²ω² − l(l+1) − α² r²),  g = (r²+a²)ω − am,

and the channel-dependent part `Δ (−l(l+1) + a²(ω²−α²) C)` of the last two terms forms one term.
"""
function _horizon_terms(background::KerrBackground, mass_coupling, frequency, m, ls, mixing, p, order::Int)
    a, gap, rplus = background.a, background.horizon_gap, background.rplus
    g_poly = [2rplus * frequency - a * m, 2rplus * frequency, frequency]
    delta_poly = [0, gap, 1]
    common_constant = 2a * m * frequency - a^2 * frequency^2
    common = _polyadd(_polymul(g_poly, g_poly),
                      _polymul(delta_poly, [common_constant - mass_coupling^2 * rplus^2, -2 * mass_coupling^2 * rplus, -mass_coupling^2]))
    orders = 0:order
    T = promote_type(typeof(p), eltype(common), eltype(mixing))
    identity_map = Matrix{T}(I, length(ls), length(ls))
    channel_map = _channel_map(T, [-l * (l + 1) for l in ls], a^2 * (frequency^2 - mass_coupling^2), mixing)
    return [SeriesTerm{T}([gap^2, 2gap, 1], 0, (p .+ orders) .* (p .+ orders .- 1), identity_map),
            SeriesTerm{T}([gap^2, 3gap, 2], 0, p .+ orders, identity_map),
            SeriesTerm{T}(common, 0, ones(T, order + 1), identity_map),
            SeriesTerm{T}(delta_poly, 0, ones(T, order + 1), channel_map)]
end

"""
    endpoint_values(series::HorizonSeries, background, s)

Profile `R` and tortoise slope `Q = dR/dr*` of each channel at offset `s`,
from `R = s^p P(s)` and `dR/dr* = (Δ/(r²+a²)) s^p (p P/s + P')`.
"""
function endpoint_values(series::HorizonSeries, background::KerrBackground, s)
    coefficients = series.coefficients
    orders = 0:size(coefficients, 1)-1
    regular_part = transpose(coefficients) * (s .^ orders)
    regular_slope = transpose(coefficients) * (orders .* s .^ (orders .- 1))
    leading_factor = exp(series.horizon_power * log(s))
    R = leading_factor .* regular_part
    Q = dr_drstar(background, s) .* leading_factor .* (series.horizon_power .* regular_part ./ s .+ regular_slope)
    return R, Q
end

# --- outer ------------------------------------------------------------------

"""
    OuterSeries

Large-radius expansion `R_{l,q} = e^{ikr} r^ν Σ_j b^{(j)}_l r^{-j}` of one
component's channels, with `k² = ω² − α²` on the decaying/outgoing branch and
`ν = −1 + i(2ω² − α²)/k`. `coefficients[j+1, channel]`.
"""
struct OuterSeries{P<:Number, T<:Number}
    wavenumber::P
    long_range_power::P
    coefficients::Matrix{T}
end

"""
Wavenumber `k_q`: `i√(α² − ω²)` (decaying) for `Re ω < α`, else `√(ω² − α²)`
(outgoing with our `e^{-iωt}` convention).
"""
outer_wavenumber(frequency, mass_coupling) = real(frequency) < mass_coupling ? im * sqrt(mass_coupling^2 - frequency^2) : sqrt(frequency^2 - mass_coupling^2)

"""Long-range power `ν_q = −1 + i(2ω² − α²)/k`."""
outer_exponent(frequency, mass_coupling, k) = -1 + im * (2 * frequency^2 - mass_coupling^2) / k

"""
    OuterSeries(background, α, ω, m, ls, mixing, leading, order)

Build the outer expansion for a component with frequency `ω`, azimuthal number
`m`, polar indices `ls`, mixing matrix `C`, and leading coefficients.
"""
function OuterSeries(background::KerrBackground, mass_coupling, frequency, m::Int, ls::AbstractVector{Int},
                     mixing::AbstractMatrix, leading::AbstractVector, order::Int)
    k = outer_wavenumber(frequency, mass_coupling)
    long_range_power = outer_exponent(frequency, mass_coupling, k)
    terms = _outer_terms(background, mass_coupling, frequency, m, ls, mixing, k, long_range_power, order)
    OuterSeries(k, long_range_power, _solve_series(terms, leading, order, 1))
end

"""
`SeriesTerm`s of `(Δ²R'' + ΔΔ'R' + ΔB R + Δ a²(ω²−α²) C R)/r⁴` applied to `R = e^{ikr} r^ν S(u)`, `u = 1/r`,
as coefficients of `u^J`. With `δ = Δ/r²`:

    δ²[(ik+νu)²S + 2(ik+νu)S' + S'' − νu²S] + δ(2−2u)u[(ik+νu)S + S'] + β S + δu² a²(ω²−α²) C S

where `S' = dS/dr = −Σ j b_j u^{j+1}`, `S'' = Σ j(j+1) b_j u^{j+2}`, and
`β = ΔB/r⁴ = (ω(1 + a²u²) − a m u²)² + δ (c₀ u² − α²)` with `c₀ = 2amω − a²ω² − l(l+1)`. The terms in `S`, `S'`
and `S''` are grouped by derivative, and the channel-dependent `−δu² l(l+1)` joins the angular coupling.
"""
function _outer_terms(background::KerrBackground, mass_coupling, frequency, m, ls, mixing, k, long_range_power, order::Int)
    a = background.a
    delta_poly = [1, -2, a^2]                                                       # δ = Δ/r²
    delta_squared = _polymul(delta_poly, delta_poly)
    wave_poly = [im * k, long_range_power]                                          # ik + νu
    edge_poly = _polymul(delta_poly, [0, 2, -2])                                    # δ (2 − 2u) u
    g_poly = [frequency, 0, a^2 * frequency - a * m]
    common_constant = 2a * m * frequency - a^2 * frequency^2
    shifted_delta = vcat([0, 0], delta_poly)                                        # δ u²
    value_poly = _polyadd(_polyadd(_polymul(delta_squared, _polymul(wave_poly, wave_poly)), -long_range_power .* vcat([0, 0], delta_squared)),
                          _polyadd(_polymul(edge_poly, wave_poly),
                                   _polyadd(_polymul(g_poly, g_poly), _polymul(delta_poly, [-mass_coupling^2, 0, common_constant]))))
    slope_poly = _polyadd(2 .* _polymul(delta_squared, wave_poly), edge_poly)
    orders = 0:order
    T = promote_type(eltype(value_poly), eltype(mixing))
    identity_map = Matrix{T}(I, length(ls), length(ls))
    channel_map = _channel_map(T, [-l * (l + 1) for l in ls], a^2 * (frequency^2 - mass_coupling^2), mixing)
    return [SeriesTerm{T}(value_poly, 0, ones(T, order + 1), identity_map),
            SeriesTerm{T}(slope_poly, 1, -orders, identity_map),                     # S' = −Σ j b_j u^{j+1}
            SeriesTerm{T}(delta_squared, 2, orders .* (orders .+ 1), identity_map),   # S'' = Σ j(j+1) b_j u^{j+2}
            SeriesTerm{T}(shifted_delta, 0, ones(T, order + 1), channel_map)]
end

"""
    endpoint_values(series::OuterSeries, background, s)

Profile `R` and tortoise slope `Q = dR/dr*` of each channel at offset `s`, from
`R = e^{ikr} r^ν S` and `dR/dr = e^{ikr} r^ν [(ik + ν/r) S + S']`.
"""
function endpoint_values(series::OuterSeries, background::KerrBackground, s)
    r = radius(background, s)
    inverse_r = 1 / r
    coefficients = series.coefficients
    orders = 0:size(coefficients, 1)-1
    regular_part = transpose(coefficients) * (inverse_r .^ orders)
    regular_slope_dr = -transpose(coefficients) * (orders .* inverse_r .^ (orders .+ 1))
    leading_factor = exp(im * series.wavenumber * r + series.long_range_power * log(r))
    R = leading_factor .* regular_part
    Q = dr_drstar(background, s) .* leading_factor .*
        ((im * series.wavenumber + series.long_range_power * inverse_r) .* regular_part .+ regular_slope_dr)
    return R, Q
end

# --- integrated outer data for propagating components -------------------------

"""
    OuterBoundary(; data=:integrated, include_bound=false, convergence_margin=40.0, far_order=40,
                  steps_per_wavelength=48)

How the outer endpoint data of *propagating* components (`Re ω_q ≥ α`) are
produced. Bound components use the truncated `OuterSeries` unless
`include_bound` is set, in which case they are integrated the same way (inward
integration of a decaying solution is the stable direction, and the envelope
`T` stays `O(1)`, so this is safe; the wavelength in `steps_per_wavelength`
then means `2π/κ`).

* `data = :series`: the truncated series at the mesh endpoint (the original
  behaviour). The series is asymptotic, and at a near-threshold wavenumber it
  is simply inaccurate at any usable endpoint (1.9e-2 in `R` at `r = 400` for
  the `(−1, 2)` component at `a = 0.99`, `α = 0.3`; see
  `validation/outer_integration_check.jl`).
* `data = :integrated`: the series is evaluated at a far radius where it
  converges, `2 k r_far ≥ |η|² + convergence_margin` with the Coulomb parameter
  `η = −(2ω² − α²)/k`, and the *linear* outer system is integrated inward to the
  endpoint. The integration is done for the slowly varying envelope
  `T = R e^{−ikr} r^{−ν}` with a 3-stage Gauss–Legendre collocation step
  (order 6, A-stable), `steps_per_wavelength` steps per `2π/k`; with the
  default settings the endpoint data are accurate to ≈ 1e-11 relative.
  `far_order` is the truncation order of the far series.
"""
struct OuterBoundary
    data::Symbol
    include_bound::Bool
    convergence_margin::Float64
    far_order::Int
    steps_per_wavelength::Int
    function OuterBoundary(; data::Symbol=:integrated, include_bound::Bool=false, convergence_margin=40.0, far_order::Int=40,
                           steps_per_wavelength::Int=48)
        data in (:integrated, :series) || throw(ArgumentError("outer data must be :integrated or :series, got $data"))
        new(data, include_bound, Float64(convergence_margin), far_order, steps_per_wavelength)
    end
end

"""Are this component's endpoint data produced by inward integration under `settings`?"""
uses_integrated_data(settings::OuterBoundary, frequency, mass_coupling) =
    settings.data == :integrated && (settings.include_bound || is_propagating(frequency, mass_coupling))

"""Does a component with this frequency propagate to infinity (outgoing branch of `outer_wavenumber`)?"""
is_propagating(frequency, mass_coupling) = real(frequency) >= mass_coupling

"""Plain floating-point value of a possibly dual real or complex number (used only to size the integration)."""
_plain_value(x::Real) = ForwardDiff.value(x)
_plain_value(z::Complex) = complex(_plain_value(real(z)), _plain_value(imag(z)))

"""
    IntegratedOuterData

Endpoint profile `R` and tortoise slope `Q` of one component's channels obtained
by integrating the linear outer system inward from `far_radius`, where the
outer series was evaluated (`far_series_tail` is the relative size of the last
retained series term there, a convergence diagnostic), with `n_steps` Gauss
steps.
"""
struct IntegratedOuterData{T<:Number}
    profiles::Vector{T}
    slopes::Vector{T}
    far_radius::Float64
    far_series_tail::Float64
    n_steps::Int
end

# 3-stage Gauss–Legendre collocation (order 6): nodes, coefficients and weights.
const _GAUSS_NODES = (0.5 - sqrt(15.0) / 10, 0.5, 0.5 + sqrt(15.0) / 10)
const _GAUSS_COEFFICIENTS = [5/36  2/9-sqrt(15.0)/15  5/36-sqrt(15.0)/30;
                             5/36+sqrt(15.0)/24  2/9  5/36-sqrt(15.0)/24;
                             5/36+sqrt(15.0)/30  2/9+sqrt(15.0)/15  5/36]
const _GAUSS_WEIGHTS = (5/18, 4/9, 5/18)

"""
    gauss_legendre_integrate(system_matrix, start, stop, state, n_steps)

Integrate the linear system `y' = M(x) y` from `start` to `stop` in `n_steps`
equal steps of the 3-stage Gauss–Legendre method (order 6, A-stable; each step
solves one dense linear system of three times the state size). `system_matrix(x)`
returns `M(x)`. Generic in the number types, so dual numbers flow through.
"""
function gauss_legendre_integrate(system_matrix, start, stop, state::AbstractVector, n_steps::Int)
    T = promote_type(eltype(state), eltype(system_matrix(start)), typeof(start), typeof(stop))
    return _gauss_legendre_integrate(T, (matrix, x) -> (matrix .= system_matrix(x)), start, stop, state, n_steps)
end

"""
    _gauss_legendre_integrate(T, system_matrix!, start, stop, state, n_steps)

`gauss_legendre_integrate` with `system_matrix!(matrix, x)` writing `M(x)` into a matrix of element type `T`; all
work arrays are allocated once per call.
"""
function _gauss_legendre_integrate(T::Type, system_matrix!, start, stop, state::AbstractVector, n_steps::Int)
    n_steps == 0 && return state
    width = (stop - start) / n_steps
    dimension = length(state)
    stage_block(i) = (i - 1) * dimension .+ (1:dimension)
    stage_matrices = [zeros(T, dimension, dimension) for _ in _GAUSS_NODES]
    collocation = zeros(T, 3dimension, 3dimension)
    stages, increment, current = zeros(T, 3dimension), zeros(T, dimension), Vector{T}(state)
    position = start
    for _ in 1:n_steps
        for (matrix, node) in zip(stage_matrices, _GAUSS_NODES)
            system_matrix!(matrix, position + node * width)
        end
        fill!(collocation, 0)
        for index in 1:3dimension
            collocation[index, index] = 1
        end
        for i in 1:3
            view(stages, stage_block(i)) .= current
            for j in 1:3
                view(collocation, stage_block(i), stage_block(j)) .-= (width * _GAUSS_COEFFICIENTS[i, j]) .* stage_matrices[j]
            end
        end
        ldiv!(lu!(collocation), stages)
        fill!(increment, 0)
        for i in 1:3
            mul!(increment, stage_matrices[i], view(stages, stage_block(i)), width * _GAUSS_WEIGHTS[i], true)
        end
        current .+= increment
        position += width
    end
    return current
end

"""
System matrix of the envelope form of the linear outer equation, written into
`system`. With `R = W T`, `W = e^{ikr} r^ν`, `φ' = W'/W = ik + ν/r`, the radial
equation `(ΔR')' + V R = 0` (`V = diag(B_l) + a²(ω²−α²) C`) becomes

    T'' = −(2φ' + Δ'/Δ) T' − [φ'² + φ'' + (Δ'/Δ) φ' + V/Δ] T,

whose solution of interest varies on the scale of `r` rather than `1/k`. The
state is `(T, T')` over the component's channels.
"""
function _envelope_system_matrix!(system::AbstractMatrix, background::KerrBackground, mass_coupling, frequency, m, ls, mixing, k,
                                  long_range_power, r)
    a = background.a
    n_chan = length(ls)
    horizon_function = r^2 - 2r + a^2
    log_delta_slope = (2r - 2) / horizon_function
    phase_slope = im * k + long_range_power / r
    phase_curvature = -long_range_power / r^2
    fill!(system, 0)
    for channel in 1:n_chan
        system[channel, n_chan + channel] = 1
        system[n_chan + channel, n_chan + channel] = -(2phase_slope + log_delta_slope)
    end
    scalar_part = phase_slope^2 + phase_curvature + log_delta_slope * phase_slope
    coupling_strength = a^2 * (frequency^2 - mass_coupling^2)
    for (row, l) in enumerate(ls)
        B = ((r^2 + a^2) * frequency - a * m)^2 / horizon_function + 2a * m * frequency - a^2 * frequency^2 -
            mass_coupling^2 * r^2 - l * (l + 1)
        system[n_chan + row, row] -= scalar_part + B / horizon_function
        for column in 1:n_chan
            system[n_chan + row, column] -= coupling_strength * mixing[row, column] / horizon_function
        end
    end
    return system
end

"""
    integrated_outer_values(background, α, ω, m, ls, mixing, leading, settings, endpoint_offset)

Endpoint data of a propagating component by the far-series-plus-inward-
integration scheme described under `OuterBoundary`. The far radius and step
count are chosen from plain values, so dual numbers in the frequency, leading
coefficients or background only flow through the arithmetic.
"""
function integrated_outer_values(background::KerrBackground, mass_coupling, frequency, m::Int, ls::AbstractVector{Int},
                                 mixing::AbstractMatrix, leading::AbstractVector, settings::OuterBoundary, endpoint_offset)
    k = outer_wavenumber(frequency, mass_coupling)
    long_range_power = outer_exponent(frequency, mass_coupling, k)
    coulomb_parameter = -(2frequency^2 - mass_coupling^2) / k
    endpoint_radius = radius(background, endpoint_offset)

    plain_k = abs(_plain_value(k))      # |k| = κ for a bound component
    plain_endpoint = _plain_value(endpoint_radius)
    far_radius = max(plain_endpoint, (abs2(_plain_value(coulomb_parameter)) + settings.convergence_margin) / (2plain_k))
    n_steps = ceil(Int, (far_radius - plain_endpoint) * plain_k * settings.steps_per_wavelength / (2pi))
    n_steps == 0 && (far_radius = plain_endpoint)

    far_series = OuterSeries(background, mass_coupling, frequency, m, ls, mixing, leading, settings.far_order)
    far_offset = n_steps == 0 ? endpoint_offset : far_radius - background.rplus
    far_profiles, far_slopes = endpoint_values(far_series, background, far_offset)
    n_steps == 0 && return IntegratedOuterData(far_profiles, far_slopes, far_radius, _series_tail(far_series, far_radius), 0)

    far_r = radius(background, far_offset)
    envelope_factor(r) = exp(im * k * r + long_range_power * log(r))
    phase_slope(r) = im * k + long_range_power / r
    envelope = far_profiles ./ envelope_factor(far_r)
    envelope_slope = far_slopes ./ (dr_drstar(background, far_offset) * envelope_factor(far_r)) .- phase_slope(far_r) .* envelope
    system_matrix!(system, r) = _envelope_system_matrix!(system, background, mass_coupling, frequency, m, ls, mixing, k,
                                                         long_range_power, r)
    initial_state = vcat(envelope, envelope_slope)
    T = promote_type(eltype(initial_state), typeof(phase_slope(far_r)), typeof(frequency), typeof(endpoint_radius), eltype(mixing))
    state = _gauss_legendre_integrate(T, system_matrix!, far_r, endpoint_radius, initial_state, n_steps)

    n_chan = length(ls)
    profiles = envelope_factor(endpoint_radius) .* state[1:n_chan]
    slopes = dr_drstar(background, endpoint_offset) .* envelope_factor(endpoint_radius) .*
             (state[n_chan + 1:end] .+ phase_slope(endpoint_radius) .* state[1:n_chan])
    return IntegratedOuterData(profiles, slopes, far_radius, _series_tail(far_series, far_radius), n_steps)
end

"""Relative size of the last retained term of the outer series at radius `r` (largest over channels)."""
function _series_tail(series::OuterSeries, r)
    coefficients = _plain_value.(series.coefficients)
    order = size(coefficients, 1) - 1
    plain_r = _plain_value(r)
    powers = plain_r .^ (-(0:order))
    return maximum(abs(coefficients[end, channel] * powers[end]) / abs(sum(coefficients[:, channel] .* powers))
                   for channel in axes(coefficients, 2))
end
