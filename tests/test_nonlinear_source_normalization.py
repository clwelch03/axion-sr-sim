"""Check the proposed source normalization against the existing Tower RHS.

Run: PYTHONPATH=src python -m unittest discover -s tests -p test_nonlinear_source_normalization.py -v
Use F=Phi instead of Phi-sin(Phi): the added source should cancel the mass
terms exactly. No evolution and no changes to the production solver.
"""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from dolan.solver import Tower
from dolan.nonlinear_source import nonlinear_source, project_real_source


class NonlinearSourceNormalizationTests(unittest.TestCase):
    def test_linear_replacement_cancels_tower_mass_terms(self):
        spin, alpha = 0.9, 0.4  # M=1: a=spin, mu=alpha.
        tower_pairs = {}
        for m in (1, 3, 5):
            parameters = dict(spin=spin, m=m, l_max=5, x_left=0, x_right=4,
                              grid_spacing_n=2, x_pml_center=-100, pml_width=10)
            tower_pairs[m] = (Tower(alpha=alpha, **parameters),
                              Tower(alpha=0, **parameters))
        radial_index = 4  # Interior node; compare only this one radius.
        reference_tower = tower_pairs[1][0]
        radius = reference_tower.r[radial_index]
        delta = reference_tower.Delta[radial_index]
        modes = [(l, m) for m, (tower, _) in tower_pairs.items()
                 for l in tower.l_values]

        cosine_theta, polar_weights = leggauss(8)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(16) / 16)[None, :]
        weights = polar_weights[:, None] * (2 * np.pi / 16)
        harmonics = np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes])
        rho_squared = radius**2 + spin**2 * cosine_theta[:, None]**2

        # Individual channels expose off-diagonal mass coupling into empty
        # neighbours; a complex mixture also checks phase and m bookkeeping.
        coefficient_cases = list(np.eye(len(modes), dtype=complex))
        coefficient_cases.append(np.array([
            1 + 0.3j, -0.2 + 0.4j, 0.1 - 0.7j,
            0.5 + 0.2j, -0.3 - 0.1j, 0.2 + 0.6j]))
        largest_relative_residual = 0.0
        for coefficients in coefficient_cases:
            real_field = 2 / radius * np.real(np.sum(
                coefficients[:, None, None] * harmonics, axis=0))
            # Proposed b_NL, but with the deliberate diagnostic substitution F=Phi.
            projected_source = project_real_source(
                real_field, radius=radius, spin=spin, alpha=alpha, delta=delta,
                harmonics=harmonics, cosine_theta=cosine_theta,
                solid_angle_weights=weights)

            # Compare the complete reusable function against the explicit
            # full-sine quadrature we derived before extracting the helper.
            expected_nonlinear = alpha**2 * radius * delta * np.sum(
                harmonics.conj() * rho_squared * (real_field - np.sin(real_field))
                * weights, axis=(1, 2))
            actual_nonlinear = nonlinear_source(
                coefficients, radius=radius, spin=spin, alpha=alpha, delta=delta,
                harmonics=harmonics, cosine_theta=cosine_theta,
                solid_angle_weights=weights)
            np.testing.assert_allclose(actual_nonlinear, expected_nonlinear,
                                       atol=1e-14, rtol=1e-12)

            mass_contributions = []
            for m, (massive, massless) in tower_pairs.items():
                channel_indices = [i for i, (_, mode_m) in enumerate(modes) if mode_m == m]
                psi = np.repeat(coefficients[channel_indices, None], massive.num_points, axis=1)
                state = (psi, np.zeros_like(psi), np.zeros_like(psi))
                # Changing only mu isolates the mass contribution to pi_dot.
                acceleration_difference = (massive.rhs(state)[1]
                                           - massless.rhs(state)[1])[:, radial_index]
                acceleration_matrix = (massive.Sigma2_0[radial_index]
                                       * np.eye(massive.n_channels)
                                       + massive.a2Delta[radial_index] * massive.c_matrix)
                # Eq. (25) needs b, not pi_dot: undo the common A^{-1} solve.
                mass_contributions.extend(acceleration_matrix @ acceleration_difference)

            residual = np.asarray(mass_contributions) + projected_source
            relative_residual = np.max(np.abs(residual)) / np.max(np.abs(projected_source))
            self.assertLess(relative_residual, 2e-13)
            # Absolute check also covers initially empty neighbour channels.
            np.testing.assert_allclose(residual, 0, atol=1e-12, rtol=0)
            largest_relative_residual = max(largest_relative_residual, relative_residual)
        print(f"\nMass cancellation at r={radius:.6f}: seven coefficient cases")
        print(f"Largest relative residual: {largest_relative_residual:.3e}")


if __name__ == "__main__":
    unittest.main()
