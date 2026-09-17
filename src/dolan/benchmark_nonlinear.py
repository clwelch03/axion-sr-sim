"""Short CPU timing/profiling probe; extrapolations exclude output and setup.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 for a single-thread baseline.
Use the same arguments on a compute node before budgeting a production run.
"""

import argparse
import cProfile
import os
import pstats
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from dolan.nonlinear_source import nonlinear_source_on_grid
from dolan.nonlinear_tower import coupled_rk4_step
from dolan.solver import Tower


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--l-max', type=int, default=5)
    parser.add_argument('--polar-count', type=int, default=32)
    parser.add_argument('--azimuthal-count', type=int, default=64)
    parser.add_argument('--grid-spacing-n', type=int, default=6)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()
    if any(value <= 0 for value in (args.l_max, args.polar_count,
                                    args.azimuthal_count, args.grid_spacing_n, args.repeats)):
        parser.error('counts and grid-spacing-n must be positive')
    towers = {m: Tower(spin=.99, alpha=.42, m=m, l_max=args.l_max,
                       x_left=-200, x_right=511, grid_spacing_n=args.grid_spacing_n,
                       x_pml_center=-100, pml_width=10)
              for m in range(1, args.l_max + 1, 2)}
    reference = towers[1]
    modes = [(l, m) for m, tower in towers.items() for l in tower.l_values]
    cosine_theta, weights = leggauss(args.polar_count)
    theta = np.arccos(cosine_theta)[:, None]
    azimuth = (2 * np.pi * np.arange(args.azimuthal_count) / args.azimuthal_count)[None, :]
    angular = dict(
        harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes]),
        cosine_theta=cosine_theta,
        solid_angle_weights=weights[:, None] * (2 * np.pi / args.azimuthal_count))
    states = {m: tuple(np.zeros((tower.n_channels, tower.num_points), complex)
                       for _ in range(3)) for m, tower in towers.items()}
    profile = np.exp(-((reference.x - 30) / 10)**2)
    profile[:2] = profile[-2:] = 0
    peak = np.max(2 * np.sqrt(3 / (8 * np.pi)) * profile / reference.r)
    states[1][0][0] = profile / peak
    states[1][1][0] = -.4j * states[1][0][0]
    coefficients = np.concatenate([state[0] for state in states.values()])
    dt = .8 * reference.dx
    step_count = int(np.ceil(1e5 / dt))
    print(f'{reference.num_points} radial nodes, {len(modes)} channels, '
          f'{args.polar_count}x{args.azimuthal_count} angular points; dt={dt:g}M')
    print('Thread environment:', {name: os.environ.get(name, 'unset') for name in
                                  ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS')})

    def nonlinear_step():
        return coupled_rk4_step(towers, states, dt, **angular)

    operations = {
        'source': lambda: nonlinear_source_on_grid(
            coefficients, radii=reference.r, spin=reference.spin, alpha=reference.alpha,
            delta=reference.Delta, **angular),
        'linear RK4': lambda: coupled_rk4_step(
            towers, states, dt, include_nonlinearity=False, **angular),
        'nonlinear RK4': nonlinear_step,
    }
    for label, operation in operations.items():
        operation()  # Warm caches and libraries before timing.
        elapsed = []
        for _ in range(args.repeats):
            start = perf_counter()
            operation()
            elapsed.append(perf_counter() - start)
        seconds = float(np.median(elapsed))
        print(f'{label}: median={seconds:.6f}s, range=[{min(elapsed):.6f}, {max(elapsed):.6f}]s')
        if 'RK4' in label:
            print(f'  {step_count} steps to 1e5M: {seconds * step_count / 3600:.2f} hours '
                  '(fixed-state timing estimate; no I/O)')
    if args.profile:
        profiler = cProfile.Profile()
        profiler.runcall(nonlinear_step)
        pstats.Stats(profiler).strip_dirs().sort_stats('cumulative').print_stats(15)


if __name__ == '__main__':
    main()
