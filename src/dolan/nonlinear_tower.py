"""Combined RHS and synchronized RK4 step for real axion towers."""

import numpy as np

from dolan.nonlinear_source import nonlinear_source_on_grid


def coupled_tower_rhs(towers, states, *, harmonics, cosine_theta,
                      solid_angle_weights, include_nonlinearity=True):
    """Return {m: (psi_dot, pi_dot, chi_dot)} for {m: Tower} and {m: state}.

    All towers must share Kerr parameters, radial grid, and PML prescription.
    Only positive odd m are supported here, with the existing same-parity l
    towers. Harmonic rows MUST follow ascending m, then each tower.l_values.
    Harmonics/quadrature are precomputed by the caller, not rebuilt each RHS.
    include_nonlinearity=False bypasses the source for exact linear regression.
    """
    if not towers or set(towers) != set(states):
        raise ValueError("towers and states must have the same nonempty set of m keys")
    ordered_m = sorted(towers)
    reference = towers[ordered_m[0]]
    channel_slices = {}
    channel_count = 0
    for m in ordered_m:
        tower = towers[m]
        if m <= 0 or m % 2 != 1 or tower.m != m:
            raise ValueError("each key must match a positive odd-m tower")
        if any(getattr(tower, name) != getattr(reference, name)
               for name in ("mass", "spin", "alpha", "dx", "x_pml_center", "pml_width")):
            raise ValueError("towers must share geometry, mass coupling, and PML parameters")
        if not np.array_equal(tower.x, reference.x):
            raise ValueError("towers must share the radial grid")
        if len(states[m]) != 3 or any(
                np.shape(field) != (tower.n_channels, tower.num_points) for field in states[m]):
            raise ValueError("each state must contain three arrays matching its tower")
        channel_slices[m] = slice(channel_count, channel_count + tower.n_channels)
        channel_count += tower.n_channels

    if not include_nonlinearity:
        return {m: towers[m].rhs(states[m]) for m in ordered_m}

    coefficients = np.concatenate([states[m][0] for m in ordered_m], axis=0)
    combined_source = nonlinear_source_on_grid(
        coefficients, radii=reference.r, spin=reference.spin, alpha=reference.alpha,
        delta=reference.Delta, harmonics=harmonics, cosine_theta=cosine_theta,
        solid_angle_weights=solid_angle_weights)
    return {m: towers[m].rhs(states[m], additional_b=combined_source[channel_slices[m]])
            for m in ordered_m}


def coupled_rk4_step(towers, states, dt, *, harmonics, cosine_theta,
                     solid_angle_weights, include_nonlinearity=True):
    """Advance all towers together by dt, returning new state arrays.

    Every stage reconstructs Phi from that stage's COMPLETE multi-m state.
    This is classical RK4 for an autonomous fixed-Kerr system. The caller
    chooses a stable timestep; this function does not select a CFL limit.
    """
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and positive")

    def rhs(stage_states):
        return coupled_tower_rhs(
            towers, stage_states, harmonics=harmonics, cosine_theta=cosine_theta,
            solid_angle_weights=solid_angle_weights,
            include_nonlinearity=include_nonlinearity)

    def add(stage_derivatives, factor):
        return {m: tuple(field + factor * derivative
                         for field, derivative in zip(states[m], stage_derivatives[m]))
                for m in states}

    k1 = rhs(states)
    k2 = rhs(add(k1, dt / 2))
    k3 = rhs(add(k2, dt / 2))
    k4 = rhs(add(k3, dt))
    return {m: tuple(field + dt / 6 * (a + 2 * b + 2 * c + d)
                     for field, a, b, c, d in zip(states[m], k1[m], k2[m], k3[m], k4[m]))
            for m in states}
