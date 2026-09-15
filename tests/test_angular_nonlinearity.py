"""Second one-radius checkpoint: project Phi - sin(Phi), without evolving.

Run: PYTHONPATH=src python -m unittest discover -s tests -p test_angular_nonlinearity.py -v
The projected quantities below are q_lm = integral Y_lm^* (Phi-sin(Phi)) dOmega.
They are NOT yet the Dolan acceleration source: mu^2, radial rescaling, and
the angular metric factor must be included when deriving that source.
"""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y
import mpmath as mp

from nonlinear_source import nonlinear_source, sine_remainder


def project_dipole_source(amplitude, polar_count=12, azimuth_count=24):
    """Construct Phi=A sin(theta) cos(varphi), then project its sine remainder.

    A is the continuum peak of the dimensionless field at this radius.
    Include every nonnegative-m channel through l=5 to check forbidden modes,
    not just the channels we expect to populate.
    """
    cosine_theta, polar_weights = leggauss(polar_count)
    theta = np.arccos(cosine_theta)[:, None]
    azimuth = (2 * np.pi * np.arange(azimuth_count) / azimuth_count)[None, :]
    solid_angle_weights = polar_weights[:, None] * (2 * np.pi / azimuth_count)
    real_field = amplitude * np.sin(theta) * np.cos(azimuth)
    # These test amplitudes are large enough to use direct subtraction.
    # A production implementation will need care for cancellation near Phi=0.
    remainder = real_field - np.sin(real_field)
    return {
        (l, m): np.sum(sph_harm_y(l, m, theta, azimuth).conj()
                       * remainder * solid_angle_weights)
        for l in range(6) for m in range(l + 1)
    }


class AngularNonlinearityTests(unittest.TestCase):
    def test_reusable_source_in_small_field_limit(self):
        # Schwarzschild makes rho^2=r^2, so the analytic q_33 prediction
        # simply acquires the factor mu^2*r^3*Delta. Radius is outside r_+=2.
        radius, alpha = 7.0, 0.4
        delta = radius**2 - 2 * radius
        cosine_theta, polar_weights = leggauss(12)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(24) / 24)[None, :]
        modes = [(1, 1), (3, 1), (3, 3), (5, 1), (5, 3), (5, 5)]
        harmonics = np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes])
        geometry = dict(radius=radius, spin=0.0, alpha=alpha, delta=delta,
                        harmonics=harmonics, cosine_theta=cosine_theta,
                        solid_angle_weights=polar_weights[:, None] * (2 * np.pi / 24))
        for amplitude in (1e-4, 1e-8):
            coefficients = np.zeros(len(modes), dtype=complex)
            # This coefficient reconstructs Phi=A sin(theta) cos(varphi).
            coefficients[0] = -amplitude * radius / (2 * np.sqrt(3 / (8 * np.pi)))
            source = nonlinear_source(coefficients, **geometry)
            expected_33 = (alpha**2 * radius**3 * delta
                           * -amplitude**3 / (48 * np.sqrt(35 / (64 * np.pi))))
            np.testing.assert_allclose(source[2], expected_33, rtol=1e-9, atol=0)
        np.testing.assert_array_equal(nonlinear_source(np.zeros(len(modes)), **geometry), 0)

    def test_sine_remainder_against_high_precision(self):
        # Both signs, zero, tiny fields, and both sides of the series switch.
        fields = np.array([-1.0, -0.1, -1e-8, 0, 1e-12, 0.0999, 0.1, 0.1001, 1.0])
        with mp.workdps(80):
            expected = np.array([float(mp.mpf(float(x)) - mp.sin(mp.mpf(float(x))))
                                 for x in fields])
        np.testing.assert_allclose(sine_remainder(fields), expected, rtol=2e-13, atol=0)
        self.assertEqual(float(sine_remainder(0.0)), 0)
        with self.assertRaises(ValueError):
            sine_remainder(np.array([1 + 1j]))

    def test_cubic_coefficient_and_amplitude_scaling(self):
        amplitude = 0.02
        source = project_dipole_source(amplitude)
        half_source = project_dipole_source(amplitude / 2)
        # Y_33 = -sqrt(35/(64*pi)) sin(theta)^3 exp(3i*varphi).
        # Phi^3/6 has exp(3i*varphi) coefficient A^3 sin(theta)^3 / 48.
        expected_cubic_33 = -amplitude**3 / (48 * np.sqrt(35 / (64 * np.pi)))
        np.testing.assert_allclose(source[3, 3], expected_cubic_33,
                                   rtol=3e-5, atol=0)
        for mode in [(1, 1), (3, 3)]:
            ratio = source[mode] / half_source[mode]
            np.testing.assert_allclose(ratio, 8, rtol=3e-5, atol=0)
            print(f"\nq{mode}(A) / q{mode}(A/2) = {ratio.real:.8f}")
        print(f"q(3,3) = {source[3, 3].real:.9e}; cubic prediction = {expected_cubic_33:.9e}")

    def test_full_sine_symmetry_and_quadrature(self):
        # Finite amplitude exercises higher powers, including generation of m=5.
        source = project_dipole_source(1.0)
        refined_source = project_dipole_source(1.0, 24, 48)
        for (l, m), coefficient in source.items():
            np.testing.assert_allclose(coefficient, refined_source[l, m],
                                       atol=1e-13, rtol=0)
            if m % 2 == 0 or (l + m) % 2 != 0:
                self.assertLess(abs(coefficient), 1e-14)
        self.assertGreater(abs(source[3, 3]), 1e-3)
        self.assertGreater(abs(source[5, 5]), 1e-5)
        largest_change = max(abs(value - refined_source[mode])
                             for mode, value in source.items())
        print(f"\nLargest change on angular grid refinement: {largest_change:.3e}")


if __name__ == "__main__":
    unittest.main()
