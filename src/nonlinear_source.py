"""Axion forcing for Dolan's equation, at one radius or across a radial grid.

Returns a contribution to b in A*pi_dot=b, not an acceleration. Angular
harmonics and quadrature are supplied by the caller so they can be precomputed.
This module does not evolve fields or implement boundary conditions.
"""

import numpy as np


def sine_remainder(real_field):
    """Evaluate Phi-sin(Phi), retaining relative accuracy near Phi=0."""
    real_field = np.asarray(real_field)
    if np.iscomplexobj(real_field):
        raise ValueError("the axion potential requires a real field")
    real_field = np.asarray(real_field, dtype=float)
    small = np.abs(real_field) < 0.1
    remainder = real_field - np.sin(real_field)
    small_field = real_field[small]
    squared_field = small_field**2
    # Through Phi^11; the omitted Phi^13/13! is negligible below |Phi|=0.1.
    remainder = np.array(remainder, copy=True)
    remainder[small] = small_field**3 * (
        1 / 6 + squared_field * (-1 / 120 + squared_field * (
            1 / 5040 + squared_field * (-1 / 362880 + squared_field / 39916800))))
    return remainder


def project_real_source(real_source, *, radius, spin, alpha, delta,
                        harmonics, cosine_theta, solid_angle_weights):
    """Return mu^2*r*Delta integral Y* rho^2 real_source dOmega, in M=1 units.

    harmonics: (channel, polar point, azimuthal point), full normalized Y_lm.
    cosine_theta: (polar point,). Weights: (polar point, 1) or full angular grid,
    INCLUDING the azimuthal spacing. real_source has the angular-grid shape.
    spin=a/M, alpha=M*mu; delta is supplied from the same geometry as the tower.
    """
    if np.iscomplexobj(real_source):
        raise ValueError("source must be real on the angular grid")
    rho_squared = radius**2 + spin**2 * np.asarray(cosine_theta)[:, None]**2
    return alpha**2 * radius * delta * np.sum(
        np.asarray(harmonics).conj() * rho_squared * real_source
        * solid_angle_weights, axis=(1, 2))


def nonlinear_source(coefficients, *, radius, spin, alpha, delta,
                     harmonics, cosine_theta, solid_angle_weights):
    """Construct real Phi and return b_NL for the supplied positive-m channels.

    coefficients has shape (channel,) in the SAME order as harmonics. All
    stored m must be positive; negative-m partners are included via 2*Re.
    Coefficients describe Phi=phi/f_a, not the canonically normalized phi.
    Angular resolution must be converged for the field/amplitude of interest.
    """
    coefficients = np.asarray(coefficients, dtype=complex)
    harmonics = np.asarray(harmonics)
    if harmonics.ndim != 3 or coefficients.shape != (harmonics.shape[0],):
        raise ValueError("expected one coefficient per angular harmonic")
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be finite and positive")
    real_field = 2 / radius * np.real(np.sum(
        coefficients[:, None, None] * harmonics, axis=0))
    return project_real_source(
        sine_remainder(real_field), radius=radius, spin=spin, alpha=alpha,
        delta=delta, harmonics=harmonics, cosine_theta=cosine_theta,
        solid_angle_weights=solid_angle_weights)


def nonlinear_source_on_grid(coefficients, *, radii, spin, alpha, delta,
                             harmonics, cosine_theta, solid_angle_weights):
    """Return b_NL with shape (channel, radial node), matching coefficients.

    radii and delta are one-dimensional arrays on the shared tower grid.
    Coefficient rows and harmonic rows must share the same positive-m channel
    order. No derivatives or radial integration enter this local source.
    Process radial blocks with matrix products; temporary angular arrays are
    bounded to 256 radial slices rather than the entire simulation grid.
    Boundary freezing/PML remain the responsibility of the evolution operator.
    """
    coefficients = np.asarray(coefficients, dtype=complex)
    radii = np.asarray(radii, dtype=float)
    delta = np.asarray(delta, dtype=float)
    harmonics = np.asarray(harmonics)
    if radii.ndim != 1 or radii.size == 0 or delta.shape != radii.shape:
        raise ValueError("radii and delta must be nonempty matching 1D arrays")
    if np.any(~np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError("radii must be finite and positive")
    if np.any(~np.isfinite(delta)):
        raise ValueError("delta must be finite")
    if harmonics.ndim != 3 or coefficients.shape != (harmonics.shape[0], radii.size):
        raise ValueError("expected coefficients shaped (channel, radial node)")
    angular_shape = harmonics.shape[1:]
    flat_harmonics = harmonics.reshape(harmonics.shape[0], -1)
    weights = np.broadcast_to(solid_angle_weights, angular_shape).ravel()
    cosine_squared = np.broadcast_to(np.asarray(cosine_theta)[:, None]**2,
                                     angular_shape).ravel()
    projection = flat_harmonics.conj() * weights
    # Phi is real: two real matrix products avoid complex angular-grid arrays.
    source = np.empty_like(coefficients)
    for start in range(0, radii.size, 256):
        stop = min(start + 256, radii.size)
        block = coefficients[:, start:stop]
        radius = radii[start:stop]
        real_field = (2 / radius) * (
            flat_harmonics.real.T @ block.real - flat_harmonics.imag.T @ block.imag)
        weighted_remainder = sine_remainder(real_field) * (
            radius**2 + spin**2 * cosine_squared[:, None])
        projected = (projection.real @ weighted_remainder
                     + 1j * (projection.imag @ weighted_remainder))
        source[:, start:stop] = alpha**2 * radius * delta[start:stop] * projected
    return source
