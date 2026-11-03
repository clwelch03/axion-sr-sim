"""Independent units, conservation, and analytic-growth checks for the catalogue."""
import importlib.util
import math
from pathlib import Path
import unittest

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "runs/adiabatic/multiphase_two_cloud/scripts/estimate_qcd_populations.py"
SPEC = importlib.util.spec_from_file_location("qcd_population_estimate", PATH)
estimate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(estimate)


class PopulationEstimateChecks(unittest.TestCase):
    def test_endpoint_matches_energy_integral_and_angular_conservation(self):
        for alpha, spin, azimuthal_number in ((.2,.99,1),(.3,.9,1),(.4,.99,1),(.4,.8,2)):
            with self.subTest(alpha=alpha, spin=spin, m=azimuthal_number):
                end = estimate.saturation_endpoint(alpha, spin, azimuthal_number)
                principal_number = azimuthal_number + 1
                epsilon_from_energy = quad(lambda mass_ratio:
                    1/(alpha*(1-(alpha*mass_ratio)**2/(2*principal_number**2))),
                    end["mass_ratio"], 1, epsabs=1e-13)[0]
                self.assertAlmostEqual(epsilon_from_energy/end["epsilon"], 1, places=11)
                self.assertAlmostEqual(spin-azimuthal_number*end["epsilon"],
                                       end["spin"]*end["mass_ratio"]**2, places=12)
                frequency = end["alpha"]*(1-end["alpha"]**2/(2*principal_number**2))
                horizon_velocity = end["spin"]/(2*(1+math.sqrt(1-end["spin"]**2)))
                self.assertAlmostEqual(frequency, azimuthal_number*horizon_velocity, places=12)

    def test_weak_binding_limit_and_nonsuperradiant_case(self):
        for azimuthal_number in (1,2):
            leading = estimate.leading_saturation(1e-5,.9,azimuthal_number)
            end = estimate.saturation_endpoint(1e-5,.9,azimuthal_number)
            self.assertLess(abs(end["epsilon"]/leading-1), 1e-4)
        self.assertEqual(estimate.saturation_endpoint(.4,.8,1)["epsilon"], 0)

    def test_qcd_mass_scaling(self):
        small, large = [estimate.physical_units(.3,mass) for mass in (5,100)]
        self.assertAlmostEqual(large["occupation_scale"]/small["occupation_scale"],400)
        self.assertAlmostEqual(large["decay_constant_gev"]/small["decay_constant_gev"],20)
        self.assertAlmostEqual(large["code_energy_to_mass_fraction"]/
                               small["code_energy_to_mass_fraction"],400)
        self.assertAlmostEqual(large["axion_ev"]/small["axion_ev"],1/20)

    def test_qcd_quartic_rate_uses_coefficient_squared(self):
        qcd = estimate.transfer_rates(.3,.99,1e18,.3)
        cosine = estimate.transfer_rates(.3,.99,1e18,1)
        for qcd_rate, cosine_rate in zip(qcd,cosine):
            self.assertAlmostEqual(qcd_rate/cosine_rate,.09)

    def test_growth_exponent_matches_quadrature_and_uncoupled_limit(self):
        occupation_scale, epsilon, dominant_rate, minor_rate, pump = 1e30,.2,.7,.02,2.1
        result = estimate.minor_growth_diagnostic(epsilon,occupation_scale,dominant_rate,minor_rate,pump)
        seed_fraction = 1/occupation_scale
        integrand = lambda time: minor_rate+pump*(seed_fraction*math.exp(dominant_rate*time))**2
        log_minor = quad(integrand,0,result["dimensionless_time"],epsabs=1e-12)[0]
        self.assertAlmostEqual(result["log10_minor_occupation"],log_minor/math.log(10),places=11)
        no_transfer = estimate.minor_growth_diagnostic(epsilon,occupation_scale,dominant_rate,minor_rate,0)
        self.assertAlmostEqual(no_transfer["log10_minor_occupation"],
                               minor_rate/dominant_rate*math.log10(epsilon*occupation_scale),places=12)

    def test_damped_minor_is_not_reported_as_a_tiny_classical_prediction(self):
        damped = estimate.minor_growth_diagnostic(.1,1e30,1,-.2,1)
        self.assertEqual(damped["status"],"minor_damped_microscopic_reference")
        unstable = estimate.minor_growth_diagnostic(.1,1e30,1,-.2,100)
        self.assertEqual(unstable["status"],"minor_instability_unresolved")
        self.assertIsNone(unstable["log10_minor_over_dominant"])


if __name__ == "__main__":
    unittest.main()
