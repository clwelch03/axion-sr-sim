"""Cheap analytic checkpoints; run with PYTHONPATH=src python -m unittest discover -s tests."""

import unittest

import numpy as np

from common.hydrogenic import growth_rate


class HydrogenicRateTests(unittest.TestCase):
    def test_211_amplitude_normalization(self):
        # Fundamental l=m=1 reduces to this closed expression; the
        # occupation-number coefficient would be 1/24 instead of 1/48.
        alpha = np.array([0.01, 0.1, 0.3])
        spin = 0.99
        horizon_factor = spin - 2 * alpha * (1 + np.sqrt(1 - spin**2))
        expected = alpha**9 / 48 * horizon_factor * (
            1 - spin**2 + horizon_factor**2)
        np.testing.assert_allclose(growth_rate(alpha, spin), expected, rtol=1e-14,
                                   atol=0)

    def test_small_alpha_scaling_and_overtone(self):
        alpha = 1e-8
        self.assertAlmostEqual(float(growth_rate(alpha, 0.9) / alpha**9),
                               0.9 / 48, places=8)
        # n_principal=3 versus 2: analytic radial normalization ratio 256/729.
        self.assertAlmostEqual(float(growth_rate(alpha, 0.9, radial_overtone=1)
                                     / growth_rate(alpha, 0.9)), 256 / 729)

    def test_threshold_and_broadcasting(self):
        spin = np.array([0.5, 0.99])
        threshold = spin / (2 * (1 + np.sqrt(1 - spin**2)))
        rates = growth_rate(np.array([0.9, 1, 1.1])[:, None] * threshold, spin)
        self.assertEqual(rates.shape, (3, 2))
        self.assertTrue(np.all(rates[0] > 0))
        np.testing.assert_allclose(rates[1], 0, atol=1e-20)
        self.assertTrue(np.all(rates[2] < 0))
        self.assertLess(growth_rate(0.1, 0), 0)
        self.assertEqual(growth_rate(0, 0.9), 0)

    def test_invalid_parameters(self):
        for kwargs in ({"radial_overtone": -1}, {"l": 0, "m": 1}, {"m": 0.5}):
            with self.assertRaises(ValueError):
                growth_rate(0.1, 0.9, **kwargs)
        for alpha, spin in ((-1, 0.9), (np.nan, 0.9), (0.1, 1), (0.1, np.inf)):
            with self.assertRaises(ValueError):
                growth_rate(alpha, spin)
