#!/usr/bin/env python
"""Mode labels and linear spectral data for the multilevel scaffold.

This module deliberately contains no evolution or transition dynamics.  It
defines the channel identity that a future multimode driver will use and makes
the single-|211> reduction explicit: an active-mode set containing only
``ModeSpec(2, 1, 1)`` must reduce to the existing m=1, l=1 tower.
"""

from dataclasses import dataclass
import os
import sys


_HERE = os.path.dirname(os.path.abspath(__file__))
_LINEAR_RATES = os.path.abspath(
    os.path.join(_HERE, "..", "..", "linear_rates")
)
sys.path.insert(0, _LINEAR_RATES)

import leaver


@dataclass(frozen=True, order=True)
class ModeSpec:
    """Identity of one hydrogenic/Leaver bound channel.

    ``principal_quantum_number`` follows ``leaver.solve``: the fundamental
    211 mode has ``n=2``.  It must not be inferred from an observed radial
    overtone in a generic time-domain signal.
    """

    principal_quantum_number: int
    ell: int
    azimuthal_number: int

    def __post_init__(self):
        if self.principal_quantum_number < 1:
            raise ValueError("principal quantum number must be positive")
        if self.ell < 0:
            raise ValueError("ell must be non-negative")
        if self.principal_quantum_number <= self.ell:
            raise ValueError("principal n must satisfy n >= ell + 1")
        if abs(self.azimuthal_number) > self.ell:
            raise ValueError("must satisfy |m| <= ell")

    @property
    def label(self):
        """Human-readable label such as ``211`` or ``322``."""
        return (
            f"{self.principal_quantum_number}"
            f"{self.ell}"
            f"{self.azimuthal_number}"
        )


MODE_211 = ModeSpec(2, 1, 1)
MODE_100 = ModeSpec(1, 0, 0)
MODE_322 = ModeSpec(3, 2, 2)


def active_mode_set(*mode_specs):
    """Return a validated immutable active-channel tuple.

    A future multimode driver should accept this explicit set rather than
    silently evolving every channel known to a catalog.  In particular,
    ``active_mode_set(MODE_211)`` is the required reduction configuration for
    comparison with the existing single-mode/tower implementation.
    """
    unique_modes = tuple(dict.fromkeys(mode_specs))
    if not unique_modes:
        raise ValueError("at least one active mode is required")
    return unique_modes


def leaver_frequency(mode_spec, alpha, dimensionless_spin,
                     spheroidal=False):
    """Return the complex Leaver frequency for one catalogued mode."""
    return complex(leaver.solve(
        alpha,
        dimensionless_spin,
        l=mode_spec.ell,
        m=mode_spec.azimuthal_number,
        n=mode_spec.principal_quantum_number,
        spheroidal=spheroidal,
    ))


def quartic_m_selection_rule(output_mode, first_mode, conjugate_mode,
                             third_mode):
    """Return whether a number-conserving quartic term is azimuthally allowed.

    The projected cubic has the structure ``psi_a psi_b* psi_c`` and therefore
    requires ``m_output = m_a - m_b + m_c``.
    """
    return output_mode.azimuthal_number == (
        first_mode.azimuthal_number
        - conjugate_mode.azimuthal_number
        + third_mode.azimuthal_number
    )


def reduction_configuration():
    """Return the mandatory single-|211> multimode regression configuration."""
    return active_mode_set(MODE_211)


if __name__ == "__main__":
    print("catalogued modes:", MODE_211.label, MODE_322.label)
    print("single-mode reduction:", [mode.label for mode in reduction_configuration()])
    print("211 frequency at alpha=0.42, chi=0.99:",
          leaver_frequency(MODE_211, 0.42, 0.99))
