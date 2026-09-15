"""Check radial indexing, locality, and the physical radial normalization.

Run: PYTHONPATH=src python -m unittest discover -s tests -p test_nonlinear_source_grid.py -v
"""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from nonlinear_source import nonlinear_source, nonlinear_source_on_grid


class NonlinearSourceGridTests(unittest.TestCase):
    def setUp(self):
        cosine_theta, polar_weights = leggauss(12)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(24) / 24)[None, :]
        self.modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
        self.angular = dict(
            harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in self.modes]),
            cosine_theta=cosine_theta,
            solid_angle_weights=polar_weights[:, None] * (2 * np.pi / 24))
        self.radii = np.array([2.5, 3.0, 7.0, 20.0])

    def test_radial_slices_and_locality(self):
        random = np.random.default_rng(42)
        coefficients = (random.normal(size=(6, 4)) + 1j * random.normal(size=(6, 4)))
        geometry = dict(radii=self.radii, spin=0.9, alpha=0.4,
                        delta=self.radii**2 - 2 * self.radii + 0.9**2)
        source = nonlinear_source_on_grid(coefficients, **geometry, **self.angular)
        self.assertEqual(source.shape, coefficients.shape)
        for index, radius in enumerate(self.radii):
            reference = nonlinear_source(
                coefficients[:, index], radius=radius, spin=0.9, alpha=0.4,
                delta=geometry['delta'][index], **self.angular)
            np.testing.assert_allclose(source[:, index], reference, rtol=2e-13, atol=1e-15)
        modified = coefficients.copy()
        modified[:, 1] = 0
        modified_source = nonlinear_source_on_grid(modified, **geometry, **self.angular)
        np.testing.assert_array_equal(modified_source[:, 1], 0)
        np.testing.assert_array_equal(modified_source[:, [0, 2, 3]], source[:, [0, 2, 3]])

    def test_cubic_radial_normalization(self):
        # In Schwarzschild, fixed psi_11=c gives Phi=A(r) sin(theta) cos(varphi)
        # with A(r)=-2*c*sqrt(3/(8*pi))/r. Thus b_33 is proportional to Delta
        # at cubic order: r^3 in the geometric prefactor cancels A(r)^3.
        coefficient = 1e-5
        coefficients = np.zeros((6, 4), dtype=complex)
        coefficients[0] = coefficient
        delta = self.radii**2 - 2 * self.radii
        source = nonlinear_source_on_grid(
            coefficients, radii=self.radii, spin=0, alpha=0.4, delta=delta, **self.angular)
        expected = (0.4**2 * (2 * coefficient * np.sqrt(3 / (8 * np.pi)))**3
                    / (48 * np.sqrt(35 / (64 * np.pi))))
        np.testing.assert_allclose(source[3] / delta, expected, rtol=1e-10, atol=0)

    def test_mismatched_radial_shapes_rejected(self):
        for coefficients, delta in ((np.zeros((4, 6)), np.ones(4)),
                                    (np.zeros((6, 4)), np.ones(3))):
            with self.assertRaises(ValueError):
                nonlinear_source_on_grid(coefficients, radii=self.radii, spin=0.9,
                                         alpha=0.4, delta=delta, **self.angular)


if __name__ == "__main__":
    unittest.main()
