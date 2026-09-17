"""RHS-only checks: no RK4 steps or cloud evolutions."""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from dolan.nonlinear_source import nonlinear_source_on_grid
from dolan.nonlinear_tower import coupled_tower_rhs
from dolan.solver import Tower


class CoupledTowerRHSTests(unittest.TestCase):
    def setUp(self):
        # Include active PML nodes and an exterior region; insertion order is
        # deliberately different from the documented harmonic row order.
        self.towers = {m: Tower(spin=0.9, alpha=0.4, m=m, l_max=5,
                               x_left=-105, x_right=5, grid_spacing_n=1,
                               x_pml_center=-100, pml_width=10) for m in (5, 1, 3)}
        cosine_theta, polar_weights = leggauss(12)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(24) / 24)[None, :]
        self.modes = [(l, m) for m in sorted(self.towers) for l in self.towers[m].l_values]
        self.angular = dict(
            harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in self.modes]),
            cosine_theta=cosine_theta,
            solid_angle_weights=polar_weights[:, None] * (2 * np.pi / 24))
        random = np.random.default_rng(10)
        self.states = {m: tuple(0.3 * (random.normal(size=(tower.n_channels, tower.num_points))
                                      + 1j * random.normal(size=(tower.n_channels, tower.num_points)))
                                for _ in range(3)) for m, tower in self.towers.items()}

    def test_disabled_source_is_exactly_linear(self):
        combined = coupled_tower_rhs(self.towers, self.states,
                                     include_nonlinearity=False, **self.angular)
        for m, tower in self.towers.items():
            for actual, expected in zip(combined[m], tower.rhs(self.states[m])):
                np.testing.assert_array_equal(actual, expected)

    def test_source_acceleration_pml_and_boundaries(self):
        combined = coupled_tower_rhs(self.towers, self.states, **self.angular)
        reference = self.towers[1]
        source = nonlinear_source_on_grid(
            np.concatenate([self.states[m][0] for m in sorted(self.towers)]),
            radii=reference.r, spin=reference.spin, alpha=reference.alpha,
            delta=reference.Delta, **self.angular)
        channel_offset = 0
        for m in sorted(self.towers):
            tower = self.towers[m]
            linear = tower.rhs(self.states[m])
            np.testing.assert_array_equal(combined[m][0], linear[0])
            np.testing.assert_array_equal(combined[m][2], linear[2])
            source_slice = source[channel_offset:channel_offset + tower.n_channels]
            for node in range(2, tower.num_points - 2):
                matrix = (tower.Sigma2_0[node] * np.eye(tower.n_channels)
                          + tower.a2Delta[node] * tower.c_matrix)
                expected_increment = np.linalg.solve(matrix, source_slice[:, node])
                np.testing.assert_allclose(combined[m][1][:, node] - linear[1][:, node],
                                           expected_increment, atol=2e-14, rtol=1e-10)
            for field in combined[m]:
                np.testing.assert_array_equal(field[:, :2], 0)
                np.testing.assert_array_equal(field[:, -2:], 0)
            channel_offset += tower.n_channels

    def test_m1_drives_initially_empty_m3(self):
        states = {m: tuple(np.zeros_like(field) for field in state)
                  for m, state in self.states.items()}
        states[1][0][0, 2:-2] = 1.0
        nonlinear = coupled_tower_rhs(self.towers, states, **self.angular)
        linear = coupled_tower_rhs(self.towers, states, include_nonlinearity=False, **self.angular)
        np.testing.assert_array_equal(linear[3][1], 0)
        self.assertGreater(np.max(np.abs(nonlinear[3][1][:, 2:-2])), 1e-6)

    def test_mismatched_towers_rejected(self):
        self.towers[3].alpha = 0.5
        with self.assertRaises(ValueError):
            coupled_tower_rhs(self.towers, self.states, **self.angular)


if __name__ == "__main__":
    unittest.main()
