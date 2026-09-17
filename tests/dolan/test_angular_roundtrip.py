"""First angular checkpoint at ONE radius; no evolution or nonlinear source.

Run: python -m unittest discover -s tests -p test_angular_roundtrip.py -v
Requires SciPy >= 1.15 for sph_harm_y(l, m, polar_angle, azimuth).
Conventions: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html
"""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y


class AngularRoundTripTests(unittest.TestCase):
    def test_real_field_at_one_radius(self):
        # Positive-m channels for the equatorially symmetric, odd-m sector.
        modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
        coefficients = np.array([
            1 + 0.3j, -0.2 + 0.4j, 0.1 - 0.7j,
            0.5 + 0.2j, -0.3 - 0.1j, 0.2 + 0.6j])
        radius = 7.0

        # dOmega = d(cos theta) dvarphi. No duplicated azimuthal endpoint.
        # These counts resolve products of harmonics through l=5, |m|=5;
        # they are NOT a resolution prescription for the sine nonlinearity.
        cosine_theta, polar_weights = leggauss(8)
        theta = np.arccos(cosine_theta)
        azimuth = 2 * np.pi * np.arange(16) / 16
        solid_angle_weights = polar_weights[:, None] * (2 * np.pi / 16)
        harmonics = np.array([
            sph_harm_y(l, m, theta[:, None], azimuth[None, :])
            for l, m in modes])  # shape: (channel, polar point, azimuthal point)

        # Independent normalization/phase check: known Y_11, including its sign.
        analytic_y11 = (-np.sqrt(3 / (8 * np.pi)) * np.sin(theta[:, None])
                        * np.exp(1j * azimuth[None, :]))
        np.testing.assert_allclose(harmonics[0], analytic_y11, atol=1e-14, rtol=0)
        # With psi_11=1, the physical field must be this real dipole divided by r.
        analytic_dipole = (-np.sqrt(3 / (2 * np.pi)) / radius
                          * np.sin(theta[:, None]) * np.cos(azimuth[None, :]))
        np.testing.assert_allclose(2 * harmonics[0].real / radius,
                                   analytic_dipole, atol=1e-14, rtol=0)

        # Reconstruct Phi = (2/r) Re sum_{m>0,l} psi_lm Y_lm.
        positive_m_sum = np.sum(coefficients[:, None, None] * harmonics, axis=0)
        real_field = 2 * positive_m_sum.real / radius

        # Check that explicit negative-m reconstruction gives the same field.
        full_sum = positive_m_sum.copy()
        for coefficient, (l, m) in zip(coefficients, modes):
            negative_coefficient = (-1)**m * coefficient.conjugate()
            full_sum += negative_coefficient * sph_harm_y(
                l, -m, theta[:, None], azimuth[None, :])
        np.testing.assert_allclose(full_sum / radius, real_field, atol=1e-14, rtol=0)

        # Project: psi_lm = r integral Y_lm^* Phi dOmega (NO extra factor of 2).
        recovered = radius * np.sum(
            harmonics.conj() * real_field * solid_angle_weights, axis=(1, 2))
        np.testing.assert_allclose(recovered, coefficients, atol=1e-13, rtol=0)

        # Orthonormality also fixes the field's angular squared norm.
        measured_norm = np.sum(real_field**2 * solid_angle_weights)
        expected_norm = 2 * np.sum(np.abs(coefficients)**2) / radius**2
        np.testing.assert_allclose(measured_norm, expected_norm, atol=1e-14, rtol=0)
        print(f"\nMaximum coefficient error: {np.max(np.abs(recovered - coefficients)):.3e}")
        print(f"Angular squared-norm error: {abs(measured_norm - expected_norm):.3e}")


if __name__ == "__main__":
    unittest.main()
