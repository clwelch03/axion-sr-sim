"""
    linear_tortoise_rhs!(derivative, state, tortoise_radius, spin, alpha,
                         frequency, modes, couplings; fundamental_m=1)

Write the linear radial derivative in tortoise coordinates into `derivative` and return it.

For `N = length(modes)`, both arrays must have length `2N`: `state[1:N]`
contains `R` and `state[N+1:2N]` contains `Q = dR/dr_*`, in `modes` order.
The arrays must not alias; `state` is otherwise unchanged.

Each mode is `(l, m)` with positive `m = n*fundamental_m`. Use the positive
integer `fundamental_m` to assign frequency `n*frequency` to each harmonic.
Supply angular `couplings` keyed by `(l, other_l, m)` and units `M = 1`.

See also [`nonlinear_tortoise_rhs!`](@ref), [`integrate_radial`](@ref).
"""
function linear_tortoise_rhs!(derivative, state, tortoise_radius, spin, alpha,
                              frequency, modes, couplings; fundamental_m=1)
    count = length(modes)
    length(state) == length(derivative) == 2count || throw(DimensionMismatch("expected [R; Q]"))
    fundamental_m isa Integer && fundamental_m > 0 || throw(ArgumentError("require positive fundamental_m"))
    geometry = tortoise_geometry(tortoise_radius, spin)
    for (index, (l, m)) in enumerate(modes)
        m > 0 && m % fundamental_m == 0 || throw(ArgumentError("require m=n*fundamental_m, n>=1"))
        # All harmonics share one pattern: n=m/m0 and omega_n=n*omega0.
        # Here frequency is omega0; the coefficient helper takes omega_n.
        coefficients = tortoise_operator_coefficients(
            geometry, spin, alpha, (m ÷ fundamental_m) * frequency, l, m)
        coupled_value = zero(eltype(state))
        for (other_index, (other_l, other_m)) in enumerate(modes)
            # Linear Kerr angular mixing preserves m. The nonlinear RHS adds
            # interactions between different m blocks through nonlinear_source.
            other_m == m || continue
            coupled_value += couplings[(l, other_l, m)] * state[other_index]
        end
        # First half: dR/dr_*=Q. Second half: the transformed radial equation.
        # No nonlinear source is included in this linear validation system.
        derivative[index] = state[count + index]
        derivative[count + index] = -coefficients.slope_coefficient * state[count + index] -
            coefficients.local_coefficient * state[index] - coefficients.coupling_scale * coupled_value
    end
    return derivative
end

"""
    boundary_state(side, position, spin, alpha, frequency, modes, couplings,
                   amplitudes; fundamental_m=1, horizon_order=3, outer_order=7)

Assemble the complex boundary state `[R; Q]` in `modes` order, with `Q = dR/dr_*`.

Choose `side = :horizon` or `:outer` and a tortoise `position`. Supply one
amplitude per `(l, m)` mode. The `frequency` is fundamental; each positive
`m = n*fundamental_m` receives frequency `n*frequency`, with positive integer
`fundamental_m`. At the outer boundary, select the bound or outgoing branch
according to each harmonic's real frequency relative to `alpha`.

Use nonnegative integer `horizon_order` and `outer_order` for the respective
linear asymptotic series. Check accuracy by varying orders and endpoints.

See also [`horizon_boundary_data`](@ref), [`outer_boundary_data`](@ref).
"""
function boundary_state(
    side, position, spin, alpha, frequency, modes, couplings, amplitudes;
    fundamental_m=1, horizon_order=3, outer_order=7,
)
    side in (:horizon, :outer) ||
        throw(ArgumentError("side must be :horizon or :outer"))
    fundamental_m isa Integer && fundamental_m > 0 ||
        throw(ArgumentError("require positive fundamental_m"))

    count = length(modes)
    length(amplitudes) == count ||
        throw(DimensionMismatch("require one amplitude per mode"))

    state_type = promote_type(typeof(complex(frequency)), eltype(amplitudes))
    state = zeros(state_type, 2count)

    for m in unique([m for (l, m) in modes])
        m > 0 && m % fundamental_m == 0 ||
            throw(ArgumentError("require m=n*fundamental_m, n>=1"))

        # Gather all angular modes sharing this harmonic frequency.
        indices = findall(mode -> mode[2] == m, modes)
        l_values = [modes[index][1] for index in indices]
        mode_frequency = (m ÷ fundamental_m) * frequency

        if side == :horizon
            data = horizon_boundary_data(
                position, spin, alpha, mode_frequency,
                l_values, m, couplings, amplitudes[indices];
                order=horizon_order,
            )
        else
            kind = real(mode_frequency) < alpha ? :bound : :outgoing
            data = outer_boundary_data(
                position, spin, alpha, mode_frequency,
                l_values, m, couplings, amplitudes[indices];
                kind, order=outer_order,
            )
        end

        # Scatter each block back into the full state.
        state[indices] = data.values
        state[count .+ indices] = data.tortoise_slopes
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
- `save_profile=false`: when true, retain both integrations and return
  `(from_horizon, from_outer, residual)`, as for nonlinear reconstruction.

See also [`solve_linear_mode`](@ref), [`nonlinear_matching_residual`](@ref).
"""
function matching_residual(frequency, horizon_amplitudes, outer_amplitudes,
                           spin, alpha, l_values, m, couplings;
                           step=0.05, inner=-100.0, outer=100.0, match=5.0,
                           save_profile=false)
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
        rhs!, vcat(horizon.values, horizon.tortoise_slopes), inner, match; step, save_profile)
    from_outer = integrate_radial(
        rhs!, vcat(infinity.values, infinity.tortoise_slopes), outer, match; step, save_profile)

    # Requiring both halves to vanish joins the profiles and their slopes.
    if save_profile
        residual = from_horizon.states[:, end] - from_outer.states[:, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon - from_outer
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
- `potential=:cosine`: select `:cosine`, `:quartic`, or `:sextic`.
- `step=0.05`: positive finite maximum RK4 spacing in tortoise radius.
- `inner=-100.0`: inner endpoint in tortoise radius.
- `outer=100.0`: outer endpoint in tortoise radius.
- `match=5.0`: matching coordinate, with `inner < match < outer`.
- `horizon_order=3`: nonnegative integer order of the horizon series.
- `outer_order=7`: nonnegative integer order of the outer series.
- `save_profile=false`: when true, return `(from_horizon, from_outer, residual)`.
  Each branch contains `positions` and `states` from its boundary to `match`.
  Keep both matching-point states so a discontinuity is visible, not averaged.

See also [`solve_nonlinear_mode`](@ref), [`matching_residual`](@ref).
"""
function nonlinear_matching_residual(
    frequency, horizon_amplitudes, outer_amplitudes,
    spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    fundamental_m=1, potential=:cosine, step=0.05,
    inner=-100.0, outer=100.0, match=5.0,
    horizon_order=3, outer_order=7,
    save_profile=false,
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
        fundamental_m, potential,
    )

    from_horizon = integrate_radial(
        rhs!, horizon_state, inner, match; step, save_profile,
    )
    from_outer = integrate_radial(
        rhs!, outer_state, outer, match; step, save_profile,
    )

    # A smooth solution requires agreement in every R and every Q.
    if save_profile
        residual = from_horizon.states[:, end] - from_outer.states[:, end]
        return (; from_horizon, from_outer, residual)
    end
    return from_horizon - from_outer
end

"""
    solve_linear_mode(frequency_guess, spin, alpha, l_values, m, couplings;
                      <keyword arguments>)

Find a linear bound mode with the first outer amplitude fixed to one.

Put the normalization mode first in `l_values`; use `[1, 3, 5]` and `m = 1`
for the fundamental cloud. Solve for the frequency, all horizon amplitudes,
and the remaining outer amplitudes in units `M = 1`. Require the real frequency
to lie strictly between zero and `alpha` during Newton iteration.

Return a named tuple with `frequency`, `horizon_amplitudes`, `outer_amplitudes`,
`residual`, `residual_norm`, and `iterations`. A small mismatch at the selected
settings does not establish endpoint or step convergence, especially for the
small imaginary frequency.

# Arguments
- `horizon_guess=ones(length(l_values))`: initial horizon amplitudes in `l_values` order.
- `outer_tail_guess=zeros(length(l_values)-1)`: initial free outer amplitudes.
- `tolerance=1e-10`: positive finite matching-residual norm tolerance.
- `max_iterations=15`: nonnegative integer Newton iteration limit.
- `difference_step=1e-6`: positive finite relative scale for centered differences.
- `integration_options...`: keywords for [`matching_residual`](@ref): `step=0.05`,
  `inner=-100.0`, `outer=100.0`, and `match=5.0`.

See also [`solve_nonlinear_mode`](@ref).
"""
function solve_linear_mode(frequency_guess, spin, alpha, l_values, m, couplings;
                           horizon_guess=ones(length(l_values)),
                           outer_tail_guess=zeros(length(l_values) - 1),
                           tolerance=1e-10, max_iterations=15,
                           difference_step=1e-6, integration_options...)
    count = length(l_values)
    count > 0 && length(horizon_guess) == count && length(outer_tail_guess) == count - 1 ||
        throw(DimensionMismatch("require N horizon and N-1 free outer amplitudes"))

    initial = vcat(complex(frequency_guess), complex.(horizon_guess), complex.(outer_tail_guess))
    residual = state -> matching_residual(
        state[1], state[2:count+1], vcat(one(state[1]), state[count+2:end]),
        spin, alpha, l_values, m, couplings; integration_options...)

    # Keep Newton trials on the positive-frequency bound branch.
    admissible = state -> all(isfinite, state) && 0 < real(state[1]) < alpha
    solution = complex_newton(residual, initial;
        tolerance, max_iterations, difference_step, admissible)

    return (; frequency=solution.state[1],
              horizon_amplitudes=solution.state[2:count+1],
              outer_amplitudes=vcat(one(solution.state[1]), solution.state[count+2:end]),
              residual=solution.residual, residual_norm=norm(solution.residual),
              iterations=solution.iterations)
end

"""
    scaled_nonlinear_matching_residual(
        state, amplitude, spin, alpha, modes, couplings,
        polar_values, azimuthal_values, grid; kwargs...,
    )

Evaluate the nonlinear matching residual divided by `amplitude = A0`.

Supply the complex solver state containing the fundamental frequency,
the horizon amplitudes divided by A0, and the free outer amplitudes
divided by A0, in that order. The fundamental outer amplitude is fixed
to A0 and is not an independent state entry.

Return the complex scaled matching residual. Forward `kwargs...` to
[`nonlinear_matching_residual`](@ref); use `save_profile=false`.

See also [`solve_nonlinear_mode`](@ref).
"""
function scaled_nonlinear_matching_residual(
    state, amplitude, spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid; kwargs...,
)
    mode_count = length(modes)

    return nonlinear_matching_residual(
        state[1],
        amplitude .* state[2:mode_count+1],
        amplitude .* vcat(one(state[1]), state[mode_count+2:end]),
        spin, alpha, modes, couplings,
        polar_values, azimuthal_values, grid;
        kwargs...,
    ) ./ amplitude
end

"""
    solve_nonlinear_mode(
    frequency_guess, amplitude, spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    horizon_guess, outer_tail_guess,
    fundamental_m=1, tolerance=1e-10,
    max_iterations=15, difference_step=1e-6,
    progress_callback=nothing, stagnation_limit=8, minimum_progress=1e-3,
    jacobian_backend=:finite, integration_options...,
)

Solve nonlinear matching at the fixed positive outer amplitude `amplitude = A0`.

Put the fundamental normalization mode first in `modes`. Supply physical
boundary amplitudes as guesses; Newton scales them by `A0` and applies its
tolerance to `residual/A0`. Require positive finite real `amplitude` and keep
the fundamental frequency on the bound branch in units `M = 1`.

Return a named tuple with `frequency`, physical `horizon_amplitudes` and
`outer_amplitudes`, physical `residual` and `residual_norm`, `scaled_residual_norm`,
and `iterations`. The first outer amplitude remains `A0`. Check step, endpoint,
series-order, and angular convergence separately from the matching tolerance.

# Arguments
- `horizon_guess`: required vector of one physical amplitude per mode.
- `outer_tail_guess`: required vector of physical outer amplitudes for `modes[2:end]`.
- `fundamental_m=1`: positive integer azimuthal number of the fundamental.
- `tolerance=1e-10`: positive finite tolerance on the scaled residual norm.
- `max_iterations=15`: nonnegative integer Newton iteration limit, or `nothing`
  to correct until convergence or numerical failure.
- `difference_step=1e-6`: positive finite relative scale for centered differences.
- `progress_callback=nothing`: optional `(iteration, residual_norm)` callback;
  the reported norm is scaled by `A0`.
- `stagnation_limit=8`, `minimum_progress=1e-3`: stagnation controls forwarded
  to [`newton_solve`](@ref).
- `jacobian_backend=:finite`: select centered differences or `:forwarddiff`
  automatic differentiation for the real-packed matching system.
- `integration_options...`: keywords for [`nonlinear_matching_residual`](@ref):
  `step=0.05`, `inner=-100.0`, `outer=100.0`, `match=5.0`, `horizon_order=3`,
  and `outer_order=7`.

See also [`solve_linear_mode`](@ref), [`real_newton`](@ref).
"""
function solve_nonlinear_mode(
    frequency_guess, amplitude, spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    horizon_guess, outer_tail_guess,
    fundamental_m=1, tolerance=1e-10,
    max_iterations=15, difference_step=1e-6,
    progress_callback=nothing, stagnation_limit=8, minimum_progress=1e-3,
    jacobian_backend=:finite, integration_options...,
)
    amplitude isa Real && isfinite(amplitude) && amplitude > 0 ||
        throw(ArgumentError("require positive finite A0"))

    count = length(modes)
    count > 0 && length(horizon_guess) == count &&
        length(outer_tail_guess) == count - 1 ||
        throw(DimensionMismatch("require N horizon and N-1 outer guesses"))
    modes[1][2] == fundamental_m ||
        throw(ArgumentError("put the fundamental normalization mode first"))

    # Unknowns: frequency, N scaled horizon amplitudes, N-1 scaled outer ones.
    # The first scaled outer amplitude is fixed to 1, so its physical value is A0.
    initial = vcat(
        complex(frequency_guess),
        complex.(horizon_guess ./ amplitude),
        complex.(outer_tail_guess ./ amplitude),
    )

    residual = state -> scaled_nonlinear_matching_residual(
        state, amplitude, spin, alpha, modes, couplings,
        polar_values, azimuthal_values, grid;
        fundamental_m, integration_options...,
    )

    # The fundamental remains bound. Higher harmonics use their own
    # boundary branches, as selected by boundary_state.
    admissible = state ->
        all(isfinite, state) && 0 < real(state[1]) < alpha

    # Reconstruction involves conjugates, so use the real Jacobian.
    solution = real_newton(
        residual, initial;
        tolerance, max_iterations, difference_step, jacobian_backend,
        admissible, progress_callback, stagnation_limit, minimum_progress,
    )

    physical_residual = amplitude .* solution.residual
    return (
        frequency=solution.state[1],
        horizon_amplitudes=amplitude .* solution.state[2:count+1],
        outer_amplitudes=amplitude .* vcat(
            one(solution.state[1]), solution.state[count+2:end],
        ),
        residual=physical_residual,
        residual_norm=norm(physical_residual),
        scaled_residual_norm=norm(solution.residual),
        iterations=solution.iterations,
    )
end
