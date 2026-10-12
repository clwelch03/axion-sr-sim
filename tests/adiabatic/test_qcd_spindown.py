"""Independent balance and unit checks for the preliminary evolution driver."""
import unittest

import numpy as np

from adiabatic.qcd_spindown import (
    EPSILON, INITIAL_CLOUD_FRACTIONS, INITIAL_SPIN, REFERENCE_TIME,
    YEAR_SECONDS, coupled_rhs, physical_state,
)


class CoupledBalanceChecks(unittest.TestCase):
    def test_nonunit_mass_matches_two_conservation_equations(self):
        # Include all channels, secondary frequency response, and b != 1.
        # This catches both a missing torque mass factor and a log-energy
        # conversion error; neither is exposed by a fixed-BH integration.
        state = np.array([0.02, 0.08, np.log(100), np.log(10), 0.0])
        coefficients = dict(
            omega1=0.145, omega2=0.146, secondary_slope=-0.0001,
            primary_horizon=-3e-10, primary_outward=1e-14,
            secondary_horizon=-1e-14, secondary_outward=2e-17,
            m0_horizon=4e-7, m0_outward=3e-18,
            m3_horizon=-1e-11, m3_outward=3e-8,
        )
        derivative = coupled_rhs(state, coefficients) * REFERENCE_TIME / YEAR_SECONDS
        mass_ratio, _, _, energies = physical_state(state)
        primary_energy, secondary_energy = energies
        normalizations = np.array([
            primary_energy, secondary_energy,
            primary_energy**2 * secondary_energy,
            primary_energy * secondary_energy**2,
        ])
        names = ['primary', 'secondary', 'm0', 'm3']
        horizon_fluxes = normalizations * [coefficients[name + '_horizon'] for name in names]
        outward_fluxes = normalizations * [coefficients[name + '_outward'] for name in names]
        torque_ratios = np.array([
            1 / coefficients['omega1'], 2 / coefficients['omega2'],
            0, 3 / (2 * coefficients['omega2'] - coefficients['omega1']),
        ])
        matrix = np.array([
            [1, 1],
            [1 / coefficients['omega1'] - 2 * secondary_energy *
             coefficients['secondary_slope'] / coefficients['omega2']**2,
             2 / coefficients['omega2']],
        ])
        total_fluxes = horizon_fluxes + outward_fluxes
        expected_cloud_rates = np.linalg.solve(
            matrix, [-sum(total_fluxes), -np.dot(torque_ratios, total_fluxes)]
        )
        np.testing.assert_allclose(
            derivative[2:4] * mass_ratio * energies, expected_cloud_rates, rtol=2e-8
        )
        self.assertAlmostEqual(derivative[0] / (-EPSILON * sum(horizon_fluxes)), 1)
        self.assertAlmostEqual(
            derivative[1] / (-EPSILON * mass_ratio * np.dot(torque_ratios, horizon_fluxes)), 1
        )
        np.testing.assert_allclose(
            -derivative[0] + EPSILON * sum(expected_cloud_rates) + derivative[4],
            0, atol=1e-28,
        )

    def test_zero_flux_preserves_state(self):
        coefficients = dict(omega1=.145, omega2=.146, secondary_slope=-1e-5)
        for name in ('primary', 'secondary', 'm0', 'm3'):
            coefficients[name + '_horizon'] = 0
            coefficients[name + '_outward'] = 0
        np.testing.assert_array_equal(coupled_rhs(np.zeros(5), coefficients), np.zeros(5))

    def test_pure_primary_does_not_numerically_seed_secondary(self):
        coefficients = dict(omega1=.145, omega2=.146, secondary_slope=0)
        for name in ('primary', 'secondary', 'm0', 'm3'):
            coefficients[name + '_horizon'] = 0
            coefficients[name + '_outward'] = 0
        coefficients['primary_horizon'] = -3e-10
        derivative = coupled_rhs(np.zeros(5), coefficients)
        self.assertEqual(derivative[3], 0)
        self.assertGreater(derivative[0], 0)
        self.assertGreater(derivative[1], 0)
        self.assertGreater(derivative[2], 0)


if __name__ == '__main__':
    unittest.main()
