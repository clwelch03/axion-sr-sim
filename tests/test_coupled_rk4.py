"""Short integration checkpoints, not a physical cloud-growth simulation."""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from dolan.nonlinear_tower import coupled_rk4_step
from dolan.solver import Tower


class CoupledRK4Tests(unittest.TestCase):
    def setUp(self):
        self.towers = {m: Tower(spin=0.9, alpha=0.4, m=m, l_max=3,
                               x_left=0, x_right=16, grid_spacing_n=1,
                               x_pml_center=-100, pml_width=10) for m in (3, 1)}
        cosine_theta, polar_weights = leggauss(10)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(20) / 20)[None, :]
        modes = [(l, m) for m in sorted(self.towers) for l in self.towers[m].l_values]
        self.angular = dict(
            harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes]),
            cosine_theta=cosine_theta,
            solid_angle_weights=polar_weights[:, None] * (2 * np.pi / 20))
        self.states = {m: tuple(np.zeros((tower.n_channels, tower.num_points), complex)
                                for _ in range(3)) for m, tower in self.towers.items()}
        profile = 3 * np.exp(-((self.towers[1].x - 8) / 3)**2)
        profile[:2] = profile[-2:] = 0
        self.states[1][0][0] = profile
        self.states[1][1][0] = -0.3j * profile

    def test_disabled_source_matches_independent_linear_steps(self):
        # Populate all fields/towers so this regression also exercises chi.
        random = np.random.default_rng(23)
        states = {m: tuple(random.normal(size=field.shape) + 1j * random.normal(size=field.shape)
                           for field in state) for m, state in self.states.items()}
        actual = coupled_rk4_step(self.towers, states, 0.05,
                                  include_nonlinearity=False, **self.angular)
        for m, tower in self.towers.items():
            for combined, independent in zip(actual[m], tower.rk4_step(states[m], 0.05)):
                np.testing.assert_array_equal(combined, independent)

    def test_one_step_transfers_to_m3_and_preserves_input_and_boundaries(self):
        saved = {m: tuple(field.copy() for field in state) for m, state in self.states.items()}
        advanced = coupled_rk4_step(self.towers, self.states, 0.1, **self.angular)
        # Both pi_3 AND psi_3 respond within one step: later stages use the
        # velocities created at earlier stages.
        self.assertGreater(np.max(np.abs(advanced[3][1])), 1e-9)
        self.assertGreater(np.max(np.abs(advanced[3][0])), 1e-10)
        for m in self.states:
            for original, snapshot, updated in zip(self.states[m], saved[m], advanced[m]):
                np.testing.assert_array_equal(original, snapshot)
                self.assertFalse(np.shares_memory(original, updated))
                np.testing.assert_array_equal(updated[:, :2], 0)
                np.testing.assert_array_equal(updated[:, -2:], 0)

    def test_fourth_order_timestep_convergence(self):
        solutions = []
        for step_count in (4, 8, 16):
            states = self.states
            for _ in range(step_count):
                states = coupled_rk4_step(self.towers, states, 0.4 / step_count, **self.angular)
            # Measure the nonlinearly generated m=3 sector specifically.
            solutions.append(np.concatenate([field.ravel() for field in states[3]]))
        coarse_difference = np.linalg.norm(solutions[0] - solutions[1])
        fine_difference = np.linalg.norm(solutions[1] - solutions[2])
        ratio = coarse_difference / fine_difference
        self.assertGreater(ratio, 14)
        self.assertLess(ratio, 18)
        print(f"\nNonlinear m=3 timestep refinement ratio: {ratio:.5f} (RK4 expects ~16)")

    def test_small_amplitude_recovers_linear_evolution(self):
        # Initially only psi_11 is populated and is real. The exact angular
        # maximum of |Y_11| is sqrt(3/(8*pi)), so normalize the peak |Phi|
        # over our radial nodes to one without depending on angular sampling.
        peak_field = np.max(2 * np.sqrt(3 / (8 * np.pi))
                            * np.abs(self.states[1][0][0]) / self.towers[1].r)
        seed = {m: tuple(field / peak_field for field in state)
                for m, state in self.states.items()}
        amplitudes = np.array([0.2, 0.1, 0.05, 0.025])
        timestep_results = []
        for dt in (0.05, 0.025):
            step_count = round(2.0 / dt)
            # Independent old linear steppers; homogeneity lets us evolve the
            # unit-amplitude seed once and compare every nonlinear run / A0.
            linear = seed
            for _ in range(step_count):
                linear = {m: tower.rk4_step(linear[m], dt)
                          for m, tower in self.towers.items()}
            linear_psi = np.concatenate([linear[m][0].ravel() for m in sorted(linear)])
            reference_norm = np.linalg.norm(linear_psi)
            errors, generated_m3 = [], []
            for amplitude in amplitudes:
                nonlinear = {m: tuple(amplitude * field for field in state)
                             for m, state in seed.items()}
                for _ in range(step_count):
                    nonlinear = coupled_rk4_step(
                        self.towers, nonlinear, dt, **self.angular)
                normalized_psi = np.concatenate([
                    nonlinear[m][0].ravel() / amplitude for m in sorted(nonlinear)])
                errors.append(np.linalg.norm(normalized_psi - linear_psi) / reference_norm)
                generated_m3.append(np.linalg.norm(nonlinear[3][0])
                                    / (amplitude * reference_norm))
            errors, generated_m3 = np.array(errors), np.array(generated_m3)
            # Cubic forcing yields O(A0^3) field differences, or O(A0^2)
            # after normalization. Test the generated tower independently too.
            for measurements in (errors, generated_m3):
                self.assertTrue(np.all(measurements > 0))
                np.testing.assert_allclose(measurements[:-1] / measurements[1:],
                                           4, rtol=0.01, atol=0)
            timestep_results.append((errors, generated_m3))
        # Verify that the observed amplitude trend survives timestep refinement.
        np.testing.assert_allclose(timestep_results[0], timestep_results[1],
                                   rtol=1e-4, atol=0)
        print("\nSmall-amplitude recovery at t=2M, dt=0.025M:")
        print("A0       normalized field discrepancy    normalized generated m=3")
        for amplitude, error, generated in zip(amplitudes, *timestep_results[1]):
            print(f"{amplitude:<8.3f} {error:.8e}                 {generated:.8e}")
        print("Discrepancy ratios:", timestep_results[1][0][:-1] / timestep_results[1][0][1:])


if __name__ == "__main__":
    unittest.main()
