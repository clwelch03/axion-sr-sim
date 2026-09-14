"""Leading small-alpha scalar growth rate, in the convention M*omega_I.

Adapted from clwelch03/axion-superradiance, src/superradiance.py,
superradiant_instability_rate (local revision f85421a6d8df13f75d6d2fe586a6c20a4d9a7dd9).
That implementation follows Ng et al., arXiv:1908.02312, Eq. (1), for the
occupation-number rate Gamma. Here M*omega_I = M*Gamma/2 because a mode
evolves as exp(-i*omega*t), and its squared amplitude grows as exp(2*omega_I*t).

The horizon factor uses omega_R ~ mu, as in the source formula. This is an
asymptotic benchmark, not a prediction of the accurate finite-alpha threshold.
"""

import math
import numbers

import numpy as np


def growth_rate(alpha, spin, *, radial_overtone=0, l=1, m=1):
    """Return signed M*omega_I; alpha and spin support NumPy broadcasting.

    alpha = M*mu >= 0; spin = J/M**2 with abs(spin) < 1.
    radial_overtone >= 0 counts radial nodes; the principal quantum number
    used by leaver.solve is radial_overtone + l + 1.
    Negative rates describe decay and are deliberately retained.
    """
    for name, value in (("radial_overtone", radial_overtone), ("l", l), ("m", m)):
        if isinstance(value, bool) or not isinstance(value, numbers.Integral):
            raise ValueError(f"{name} must be an integer")
    if radial_overtone < 0 or l < 0 or abs(m) > l:
        raise ValueError("require radial_overtone >= 0, l >= 0, and abs(m) <= l")

    alpha, spin = np.broadcast_arrays(
        np.asarray(alpha, dtype=float), np.asarray(spin, dtype=float))
    if np.any(~np.isfinite(alpha)) or np.any(alpha < 0):
        raise ValueError("alpha must be finite and nonnegative")
    if np.any(~np.isfinite(spin)) or np.any(np.abs(spin) >= 1):
        raise ValueError("spin must be finite with abs(spin) < 1")

    principal_number = radial_overtone + l + 1
    horizon_radius_over_mass = 1 + np.sqrt(1 - spin**2)
    horizon_factor = m * spin - 2 * alpha * horizon_radius_over_mass

    # Log factorials avoid constructing very large combinatorial factors.
    log_prefactor = (
        (4 * l + 2) * math.log(2)
        + math.lgamma(2 * l + 2 + radial_overtone)
        - (2 * l + 4) * math.log(principal_number)
        - math.lgamma(radial_overtone + 1)
        + 2 * (math.lgamma(l + 1) - math.lgamma(2 * l + 1)
               - math.lgamma(2 * l + 2)))
    angular_product = np.ones_like(horizon_factor)
    for k in range(1, l + 1):
        angular_product *= k**2 * (1 - spin**2) + horizon_factor**2

    occupation_rate = (math.exp(log_prefactor) * alpha**(4 * l + 5)
                       * horizon_factor * angular_product)
    return occupation_rate / 2
