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

function _polyadd!(accumulator::AbstractVector, addend::AbstractVector)
    accumulator[1:length(addend)] .+= addend
    return accumulator
end

"""Polynomial product of `poly(v)` with every column of `coefficients`, padded or truncated to `n_rows`."""
function _polymul_columns(poly::AbstractVector, coefficients::AbstractMatrix, n_rows::Int)
    result = zeros(promote_type(eltype(poly), eltype(coefficients)), n_rows, size(coefficients, 2))
    for column in axes(coefficients, 2)
        product = _polymul(poly, coefficients[:, column])
        kept = 1:min(n_rows, length(product))
        result[kept, column] = product[kept]
    end
    return result
end

# --- shared order-by-order solver -------------------------------------------

"""
    _solve_series(apply_operator, leading, order, row_shift)

Fill a coefficient matrix `coefficients[j+1, channel]`, `j = 0:order`, such
that the residual polynomial `apply_operator(coefficients)` vanishes through
order `order + row_shift`. Unknown row `j` is determined by residual row
`j + row_shift`; the small matrix coupling the channels at that order is
extracted numerically from the operator.
"""
function _solve_series(apply_operator, leading::AbstractVector, order::Int, row_shift::Int)
    n_chan = length(leading)
    T = promote_type(eltype(leading), eltype(apply_operator(reshape(leading, 1, n_chan))))
    coefficients = zeros(T, order + 1, n_chan)
    coefficients[1, :] = leading
    for j in 1:order
        residual = apply_operator(coefficients)[j + 1 + row_shift, :]
        coupling = zeros(T, n_chan, n_chan)
        for channel in 1:n_chan
            unit = zeros(T, order + 1, n_chan)
            unit[j + 1, channel] = one(T)
            coupling[:, channel] = apply_operator(unit)[j + 1 + row_shift, :]
        end
        coefficients[j + 1, :] = -(coupling \ residual)
    end
    return coefficients
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
    apply_operator(coefficients) = _horizon_operator(background, mass_coupling, frequency, m, ls, mixing, p, coefficients)
    HorizonSeries(p, _solve_series(apply_operator, leading, order, 0))
end

"""
Residual coefficients of `Δ(ΔR')' + ΔB R + Δ a²(ω²−α²) C R` applied to
`R = s^p Σ_j c_j s^j`, as rows of `s^{p+J}`. With `Δ = s(s + gap)` every
coefficient is a polynomial in `s`.
"""
function _horizon_operator(background::KerrBackground, mass_coupling, frequency, m, ls, mixing, p, coefficients::AbstractMatrix)
    a, gap, rplus = background.a, background.horizon_gap, background.rplus
    n_orders, n_chan = size(coefficients)
    n_rows = n_orders + 4
    orders = 0:n_orders-1
    residual = zeros(promote_type(eltype(coefficients), typeof(p)), n_rows, n_chan)

    # Δ² R'' = (gap² s² + 2 gap s³ + s⁴) Σ (p+j)(p+j−1) c_j s^{p+j−2}
    residual .+= _polymul_columns([gap^2, 2gap, 1], (p .+ orders) .* (p .+ orders .- 1) .* coefficients, n_rows)
    # Δ Δ' R' = (gap² s + 3 gap s² + 2 s³) Σ (p+j) c_j s^{p+j−1}
    residual .+= _polymul_columns([gap^2, 3gap, 2], (p .+ orders) .* coefficients, n_rows)
    # ΔB R with ΔB = g² + Δ (2amω − a²ω² − l(l+1) − α² r²),  g = (r²+a²)ω − am
    g_poly = [2rplus * frequency - a * m, 2rplus * frequency, frequency]
    delta_poly = [0, gap, 1]
    for (channel, l) in enumerate(ls)
        constant = 2a * m * frequency - a^2 * frequency^2 - l * (l + 1)
        delta_B = _polyadd!(_polymul(g_poly, g_poly),
                            _polymul(delta_poly, [constant - mass_coupling^2 * rplus^2, -2 * mass_coupling^2 * rplus, -mass_coupling^2]))
        residual[:, channel] .+= _polymul_columns(delta_B, coefficients[:, channel:channel], n_rows)
    end
    # angular coupling Δ a²(ω²−α²) Σ_l' C_ll' R_l'
    coupling_strength = a^2 * (frequency^2 - mass_coupling^2)
    residual .+= _polymul_columns(delta_poly .* coupling_strength, coefficients * transpose(mixing), n_rows)
    return residual
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
    apply_operator(coefficients) = _outer_operator(background, mass_coupling, frequency, m, ls, mixing, k, long_range_power, coefficients)
    OuterSeries(k, long_range_power, _solve_series(apply_operator, leading, order, 1))
end

"""
Residual coefficients of `(Δ²R'' + ΔΔ'R' + ΔB R + Δ a²(ω²−α²) C R)/r⁴` applied
to `R = e^{ikr} r^ν S(u)`, `u = 1/r`, as rows of `u^J`. With `δ = Δ/r²`:

    δ²[(ik+νu)²S + 2(ik+νu)S' + S'' − νu²S] + δ(2−2u)u[(ik+νu)S + S'] + β S + δu² a²(ω²−α²) C S

where `S' = dS/dr = −Σ j b_j u^{j+1}`, `S'' = Σ j(j+1) b_j u^{j+2}`, and
`β = ΔB/r⁴ = (ω(1 + a²u²) − a m u²)² + δ (c₀ u² − α²)`.
"""
function _outer_operator(background::KerrBackground, mass_coupling, frequency, m, ls, mixing, k, long_range_power, coefficients::AbstractMatrix)
    a = background.a
    n_orders, n_chan = size(coefficients)
    n_rows = n_orders + 6
    orders = collect(0:n_orders-1)
    T = promote_type(eltype(coefficients), typeof(long_range_power))
    residual = zeros(T, n_rows, n_chan)

    series = coefficients
    series_dr = vcat(zeros(T, 1, n_chan), -orders .* coefficients)                  # −j b_j u^{j+1}
    series_d2r = vcat(zeros(T, 2, n_chan), orders .* (orders .+ 1) .* coefficients) # j(j+1) b_j u^{j+2}
    delta_poly = [1, -2, a^2]                                                       # δ = Δ/r²
    delta_squared = _polymul(delta_poly, delta_poly)
    wave_poly = [im * k, long_range_power]                                                         # ik + νu
    edge_poly = _polymul(delta_poly, [0, 2, -2])                                    # δ (2 − 2u) u

    residual .+= _polymul_columns(_polyadd!(_polymul(delta_squared, _polymul(wave_poly, wave_poly)),
                                            -long_range_power .* vcat([0, 0], delta_squared)), series, n_rows)
    residual .+= _polymul_columns(2 .* _polymul(delta_squared, wave_poly), series_dr, n_rows)
    residual .+= _polymul_columns(delta_squared, series_d2r, n_rows)
    residual .+= _polymul_columns(_polymul(edge_poly, wave_poly), series, n_rows)
    residual .+= _polymul_columns(edge_poly, series_dr, n_rows)
    g_poly = [frequency, 0, a^2 * frequency - a * m]
    for (channel, l) in enumerate(ls)
        constant = 2a * m * frequency - a^2 * frequency^2 - l * (l + 1)
        beta_poly = _polyadd!(_polymul(g_poly, g_poly), _polymul(delta_poly, [-mass_coupling^2, 0, constant]))
        residual[:, channel] .+= _polymul_columns(beta_poly, coefficients[:, channel:channel], n_rows)
    end
    coupling_strength = a^2 * (frequency^2 - mass_coupling^2)
    residual .+= _polymul_columns(vcat([0, 0], delta_poly) .* coupling_strength, coefficients * transpose(mixing), n_rows)
    return residual
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
