"""
    integrate_radial(rhs!, initial_state, start, stop; step=0.05)

Integrate from `start` to `stop` using fixed-step fourth-order Runge–Kutta.

The callback `rhs!(derivative, state, position)` must fill `derivative` without
changing `state`. Return only the final state as a complex floating-point array;
leave `initial_state` unchanged. Integrate in either direction and adjust the
spacing to land on `stop`.

For radial evolution, positions are tortoise coordinates and `state = [R; Q]`
with `Q = dR/dr_*`. Require positive finite `step`, which sets the maximum
spacing, not an error tolerance. Check accuracy by reducing it.
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

"""
    matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                      spin, alpha, l_values, m, couplings; <keyword arguments>)

Return the linear mismatch `[R_horizon-R_outer; Q_horizon-Q_outer]` at `match`.

Use one positive fundamental `m` and complex bound-state `frequency` in units
`M = 1`. Both amplitude vectors follow `l_values` order. Horizon and outer
series use orders three and seven, respectively. A root solve must fix one
nonzero amplitude to exclude the trivial zero field.

# Arguments
- `step=0.05`: positive finite maximum RK4 spacing in tortoise radius.
- `inner=-100.0`: inner endpoint in tortoise radius.
- `outer=100.0`: outer endpoint in tortoise radius.
- `match=5.0`: matching coordinate, with `inner < match < outer`.

See also [`solve_linear_mode`](@ref), [`nonlinear_matching_residual`](@ref).
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


"""
    nonlinear_matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                                spin, alpha, modes, couplings, polar_values,
                                azimuthal_values, grid; <keyword arguments>)

Return the nonlinear mismatch `[R_horizon-R_outer; Q_horizon-Q_outer]` at `match`.

Use the fundamental `frequency` in units `M = 1`; both amplitude vectors follow
`modes` order. Construct linear asymptotic boundaries with [`boundary_state`](@ref)
and integrate the coupled nonlinear equations from each endpoint. Supply angular
tables on `grid` for all modes. This evaluates the mismatch without fixing a
normalization or finding a frequency.

# Arguments
- `fundamental_m=1`: positive integer azimuthal number of the fundamental.
- `step=0.05`: positive finite maximum RK4 spacing in tortoise radius.
- `inner=-100.0`: inner endpoint in tortoise radius.
- `outer=100.0`: outer endpoint in tortoise radius.
- `match=5.0`: matching coordinate, with `inner < match < outer`.
- `horizon_order=3`: nonnegative integer order of the horizon series.
- `outer_order=7`: nonnegative integer order of the outer series.

See also [`solve_nonlinear_mode`](@ref), [`matching_residual`](@ref).
"""
function nonlinear_matching_residual(
    frequency, horizon_amplitudes, outer_amplitudes,
    spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    fundamental_m=1, step=0.05,
    inner=-100.0, outer=100.0, match=5.0,
    horizon_order=3, outer_order=7,
)
    inner < match < outer ||
        throw(ArgumentError("match must lie between endpoints"))

    # Assemble all harmonic blocks, using their respective frequencies
    # and ingoing, decaying, or outgoing boundary conditions.
    horizon_state = boundary_state(
        :horizon, inner, spin, alpha, frequency,
        modes, couplings, horizon_amplitudes;
        fundamental_m, horizon_order,
    )
    outer_state = boundary_state(
        :outer, outer, spin, alpha, frequency,
        modes, couplings, outer_amplitudes;
        fundamental_m, outer_order,
    )

    # Both integrations use the same coupled nonlinear equations.
    # Angular tables are reused at every RK4 evaluation.
    rhs! = (derivative, state, position) -> nonlinear_tortoise_rhs!(
        derivative, state, position, spin, alpha, frequency,
        modes, couplings, polar_values, azimuthal_values, grid;
        fundamental_m,
    )

    from_horizon = integrate_radial(
        rhs!, horizon_state, inner, match; step,
    )
    from_outer = integrate_radial(
        rhs!, outer_state, outer, match; step,
    )

    # A smooth solution requires agreement in every R and every Q.
    return from_horizon - from_outer
end

