"""Independent checks of the population catch-up timescale calculation."""
import math
from pathlib import Path
import sys
import unittest

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runs/adiabatic/multiphase_two_cloud/scripts"))
import check_qcd_catchup as catchup
from estimate_qcd_populations import growth_rate, threshold_spin


class CatchupChecks(unittest.TestCase):
    def test_rate_convention_and_threshold(self):
        for alpha in (0.05, 0.2, 0.4):
            for spin in (0.5, 0.8, 0.99):
                for m in (1, 2):
                    expected = 2*float(growth_rate(alpha, spin, l=m, m=m))
                    actual = catchup.occupation_rate(alpha, spin, m, binding=False)
                    self.assertAlmostEqual(actual/expected, 1.0, places=13)
            for m in (1, 2):
                omega = alpha*(1-alpha**2/(2*(m+1)**2))
                threshold = threshold_spin(omega, m)
                self.assertLess(abs(catchup.occupation_rate(alpha, threshold, m)), 1e-20)
                self.assertLess(catchup.occupation_rate(alpha, threshold-1e-4, m), 0)
                self.assertGreater(catchup.occupation_rate(alpha, threshold+1e-4, m), 0)

    def test_analytic_gain_against_independent_quadrature(self):
        half_life, horizon, pump, transition = 2.7, 0.13, 4.2, 0.08
        for duration in (0.1, half_life, 17.0):
            integrand = lambda time: (horizon + pump/(1+time/half_life)**2
                                     + transition/(1+time/half_life))
            expected = quad(integrand, 0, duration, epsabs=1e-12)[0]
            actual = catchup.growth_exponent_with_gw_depletion(
                duration, half_life, horizon, pump, transition)
            self.assertAlmostEqual(actual, expected, places=11)
        # At one GW half-life the pump contributes only half the frozen-rate gain.
        self.assertAlmostEqual(catchup.growth_exponent_with_gw_depletion(
            half_life, half_life, 0, pump), pump*half_life/2)

    def test_post_extraction_occupation_units(self):
        # Construct endpoints from an independent energy integral so this check
        # also runs in a fresh checkout without ignored campaign results.
        for alpha, dominant, mass_ratio in ((0.2, 1, 0.98), (0.3, 1, 0.99), (0.4, 2, 0.98)):
            principal = dominant + 1
            epsilon = quad(lambda mass: 1/(alpha*(1-(alpha*mass)**2/(2*principal**2))),
                           mass_ratio, 1, epsabs=1e-13)[0]
            final_alpha = alpha*mass_ratio
            final_frequency = final_alpha*(1-final_alpha**2/(2*principal**2))
            final_spin = threshold_spin(final_frequency, dominant)
            initial_spin = final_spin*mass_ratio**2 + dominant*epsilon
            row = dict(alpha=alpha, spin=initial_spin, dominant_m=dominant,
                       occupation_scale=1e70, time_unit_seconds=1e-4)
            point = catchup.endpoint_coordinates(row, "mass_changing")
            expected_number = epsilon*row["occupation_scale"]
            self.assertAlmostEqual(point["epsilon"]*point["occupation_scale"]/expected_number, 1.0)
            removed_angular_momentum = dominant*point["epsilon_initial_units"]
            final_angular_momentum = point["spin"]*point["mass_ratio"]**2
            self.assertAlmostEqual(initial_spin-removed_angular_momentum, final_angular_momentum)
            self.assertAlmostEqual(point["alpha"]/alpha, mass_ratio)
            self.assertAlmostEqual(point["time_unit_seconds"]/row["time_unit_seconds"], mass_ratio)

    def test_cosine_and_qcd_coupling_scaling(self):
        cosine = catchup.rates(0.2, 0.8, 1e18)
        qcd = catchup.rates(0.2, 0.8, 1e18, 0.3)
        for key in ("transfer_bh", "transfer_infinity", "relativistic_211", "relativistic_322"):
            self.assertAlmostEqual(qcd[key]/cosine[key], 0.09)
        for key in ("annihilation_211", "annihilation_322", "transition_322_to_211"):
            self.assertEqual(qcd[key], cosine[key])


if __name__ == "__main__":
    unittest.main()
