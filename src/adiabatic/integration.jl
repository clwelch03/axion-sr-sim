"""Integrate from start to stop in r_* using fixed-step RK4.

rhs!(derivative, state, position) fills the derivative without changing state.
For our radial system, state=[R; Q] with Q=dR/dr_*. Return only the final
state; initial_state is left unchanged. step is a positive maximum spacing,
not an error tolerance: accuracy must be checked by reducing it.
"""
function integrate_radial(rhs!, initial_state, start, stop; step=0.05)
    isfinite(step) && step > 0 ||
        throw(ArgumentError("step must be positive and finite"))

    # Adjust the spacing to land exactly on the destination
    # It becomes negative when integrating inward
    count = max(1, ceil(Int, abs(stop - start) / step))
    spacing = (stop - start) / count
    state = complex.(float.(initial_state))
    k1, k2, k3, k4 = (similar(state) for _ in 1:4)

    for index in 0:(count - 1)
        position = start + index * spacing

        # Four estimates of the derivative across this step
        rhs!(k1, state, position)
        rhs!(k2, state .+ (spacing / 2) .* k1, position + spacing / 2)
        rhs!(k3, state .+ (spacing / 2) .* k2, position + spacing / 2)
        rhs!(k4, state .+ spacing .* k3, position + spacing)

        state .+= (spacing / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end

    return state
end

"""Return [R_horizon-R_outer; Q_horizon-Q_outer] at the matching point.

This is a linear bound-state calculation for one positive fundamental m.
frequency is its complex frequency; amplitude vectors follow l_values order.
All positions are tortoise coordinates. The boundary series use their default
orders (3 near the horizon, 7 at infinity).

A later root solve must fix one nonzero amplitude, such as the outer l=1
coefficient, to exclude the trivial zero field. This function only evaluates
the mismatch; it neither enforces that normalization nor finds a frequency.
"""
function matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                           spin, alpha, l_values, m, couplings;
                           step=0.05, inner=-100.0, outer=100.0, match=5.0)
    inner < match < outer ||
        throw(ArgumentError("match must lie between endpoints"))

    # This calculation contains one fundamental-m block
    modes = [(l, m) for l in l_values]
    rhs! = (derivative, state, position) -> linear_tortoise_rhs!(
        derivative, state, position, spin, alpha, frequency, modes, couplings;
        fundamental_m=m)

    # Boundary series supply both R and Q at each endpoint
    horizon = horizon_boundary_data(
        inner, spin, alpha, frequency, l_values, m, couplings, horizon_amplitudes)
    infinity = outer_boundary_data(
        outer, spin, alpha, frequency, l_values, m, couplings, outer_amplitudes)

    from_horizon = integrate_radial(
        rhs!, vcat(horizon.values, horizon.tortoise_slopes), inner, match; step)
    from_infinity = integrate_radial(
        rhs!, vcat(infinity.values, infinity.tortoise_slopes), outer, match; step)

    # Requiring both halves to vanish joins the profiles and their slopes.
    return from_horizon - from_infinity
end
