"""Finite-amplitude angular sampling checks at fixed harmonic truncation.

These do not establish convergence in l_max, m_max, radius, or evolution time.
"""

import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from nonlinear_source import nonlinear_source
from nonlinear_tower import coupled_rk4_step
from solver import Tower


def angular_grid(modes, polar_count, azimuthal_count):
    cosine_theta, polar_weights = leggauss(polar_count)
    theta = np.arccos(cosine_theta)[:, None]
    azimuth = (2 * np.pi * np.arange(azimuthal_count) / azimuthal_count)[None, :]
    return dict(
        harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes]),
        cosine_theta=cosine_theta,
        solid_angle_weights=polar_weights[:, None] * (2 * np.pi / azimuthal_count))


class FiniteAmplitudeQuadratureTests(unittest.TestCase):
    def test_mixed_mode_source_polar_and_azimuthal_refinement(self):
        modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
        coefficients = np.array([1, .3 + .2j, -.2j, .25j, -.15, .1 + .1j])
        radius = 7.0
        reference_grid = angular_grid(modes, 64, 128)
        field = 2 / radius * np.real(np.sum(
            coefficients[:, None, None] * reference_grid['harmonics'], axis=0))
        # Fix the same coefficients across all quadratures; do not renormalize
        # against each grid's sampled peak, which would change the input.
        coefficients /= np.max(np.abs(field))
        geometry = dict(radius=radius, spin=.9, alpha=.4,
                        delta=radius**2 - 2 * radius + .9**2)
        for amplitude in (1.0, 3.0):
            reference = nonlinear_source(amplitude * coefficients,
                                         **geometry, **reference_grid)
            finer_reference = nonlinear_source(
                amplitude * coefficients, **geometry, **angular_grid(modes, 80, 160))
            np.testing.assert_allclose(reference, finer_reference, rtol=1e-11, atol=1e-12)
            for direction in ('polar', 'azimuthal'):
                errors = []
                for count in (8, 16, 32):
                    grid = angular_grid(modes, count if direction == 'polar' else 64,
                                        2 * count if direction == 'azimuthal' else 128)
                    source = nonlinear_source(amplitude * coefficients, **geometry, **grid)
                    errors.append(np.linalg.norm(source - reference) / np.linalg.norm(reference))
                self.assertGreater(errors[0], 100 * errors[1])
                self.assertLess(errors[2], 1e-10)
                print(f'\nPeak Phi~{amplitude:g}, {direction} source errors: {errors}')

    def test_short_finite_amplitude_evolution_quadrature_refinement(self):
        towers = {m: Tower(spin=.9, alpha=.4, m=m, l_max=3,
                           x_left=0, x_right=16, grid_spacing_n=1,
                           x_pml_center=-100, pml_width=10) for m in (1, 3)}
        modes = [(l, m) for m in sorted(towers) for l in towers[m].l_values]
        seed = {m: tuple(np.zeros((tower.n_channels, tower.num_points), complex)
                         for _ in range(3)) for m, tower in towers.items()}
        profile = np.exp(-((towers[1].x - 8) / 3)**2)
        profile[:2] = profile[-2:] = 0
        peak = np.max(2 * np.sqrt(3 / (8 * np.pi)) * profile / towers[1].r)
        seed[1][0][0] = 3 * profile / peak
        seed[1][1][0] = -.3j * seed[1][0][0]
        solutions = []
        for polar_count in (6, 12, 24, 32):
            angular = angular_grid(modes, polar_count, 2 * polar_count)
            states = seed
            for _ in range(20):
                states = coupled_rk4_step(towers, states, .025, **angular)
            # Include the newly generated m=3 sector separately so a large
            # unchanged m=1 component cannot hide its quadrature error.
            solutions.append([np.concatenate([states[m][index].ravel() for m in sector
                                              for index in (0, 1)])
                              for sector in ((1, 3), (3,))])
        for sector_index, label in enumerate(('all modes', 'm=3')):
            reference = solutions[-1][sector_index]
            errors = [np.linalg.norm(solution[sector_index] - reference)
                      / np.linalg.norm(reference) for solution in solutions[:-1]]
            self.assertGreater(errors[0], 100 * errors[1])
            self.assertLess(errors[1], 1e-8)
            self.assertLess(errors[2], 1e-11)
            print(f'\nPeak Phi=3, t=.5M, {label} evolution errors: {errors}')


if __name__ == '__main__':
    unittest.main()
