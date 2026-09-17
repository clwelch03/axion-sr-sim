"""Paired full-RK4 benchmark including Julia communication and array copies.

The production solver is unchanged: its remainder function is substituted
only inside this experiment, then restored. Julia stays alive across calls.
"""

import argparse
import cProfile
from collections import defaultdict
from functools import partial
import inspect
from pathlib import Path
import pstats
import sys
from time import perf_counter
import unittest
from unittest.mock import patch

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y

from dolan.julia_remainder_bridge import JuliaRemainder
from dolan import nonlinear_source
from dolan import nonlinear_tower
from dolan.nonlinear_tower import coupled_rk4_step
from dolan.solver import Tower


def profile_source_parts(towers, states, dt, angular, bridge, repeats):
    """Time existing source lines; nested bridge time is reported separately."""
    source = nonlinear_source.nonlinear_source_on_grid
    lines, first_line = inspect.getsourcelines(source)
    def line_number(prefix):
        return first_line + next(index for index, line in enumerate(lines)
                                 if line.strip().startswith(prefix))
    reconstruction = line_number('real_field =')
    weighting = line_number('weighted_remainder =')
    projection = line_number('projected =')
    assembly = line_number('source[:, start:stop] =')
    totals = defaultdict(float)
    bridge_calls = 0

    def timed_remainder(field):
        nonlocal bridge_calls
        start = perf_counter()
        result = bridge(field)
        totals['remainder bridge'] += perf_counter() - start
        bridge_calls += 1
        return result

    def trace_source(frame, event, argument):
        if event != 'call' or frame.f_code is not source.__code__:
            return None
        previous_time, previous_line = perf_counter(), frame.f_lineno
        def trace_line(frame, event, argument):
            nonlocal previous_time, previous_line
            if event in ('line', 'return'):
                now = perf_counter()
                if reconstruction <= previous_line < weighting:
                    category = 'reconstruction'
                elif weighting <= previous_line < projection:
                    category = 'weighting including bridge'
                elif projection <= previous_line < assembly:
                    category = 'projection'
                else:
                    category = 'source setup and assembly'
                totals[category] += now - previous_time
                previous_time, previous_line = now, frame.f_lineno
            return trace_line
        return trace_line

    elapsed = [[], []]
    with patch.object(nonlinear_source, 'sine_remainder', timed_remainder):
        expected = coupled_rk4_step(towers, states, dt, **angular)
        for repeat in range(repeats):
            for index in ((0, 1) if repeat % 2 == 0 else (1, 0)):
                totals.clear()
                bridge_calls = 0
                previous_trace = sys.gettrace()
                start = perf_counter()
                try:
                    if index:
                        sys.settrace(trace_source)
                    actual = coupled_rk4_step(towers, states, dt, **angular)
                finally:
                    sys.settrace(previous_trace)
                elapsed[index].append(perf_counter()-start)
                if index:
                    totals['weighting'] = totals.pop('weighting including bridge') - totals['remainder bridge']
                    print(f'Profiled step {repeat+1}: total={elapsed[index][-1]:.6f}s; '
                          + '; '.join(f'{name}={value:.6f}s' for name, value in totals.items())
                          + f'; bridge calls={bridge_calls}', flush=True)
                for m in towers:
                    for target, result in zip(expected[m], actual[m]):
                        np.testing.assert_array_equal(target, result)
    print('Instrumented/uninstrumented states agree exactly.', flush=True)
    for name, timings in zip(('Untraced', 'Line-traced'), elapsed):
        print(f'{name}: median={np.median(timings):.6f}s; '
              f'range={min(timings):.6f}..{max(timings):.6f}s', flush=True)


def compare_layouts(towers, states, dt, angular, bridge, repeats):
    """Include packing costs; compare source values and complete RK4 states."""
    source = nonlinear_source.nonlinear_source_on_grid
    reference = towers[1]
    geometry = dict(radii=reference.r, spin=reference.spin,
                    alpha=reference.alpha, delta=reference.Delta, **angular)
    # Exercise every channel and both coefficient components, not only the seed.
    generator = np.random.default_rng(420)
    shape = (angular['harmonics'].shape[0], reference.num_points)
    coefficients = (generator.normal(size=shape) + 1j*generator.normal(size=shape))
    coefficients *= reference.r[None, :] / 6
    with patch.object(nonlinear_source, 'sine_remainder', bridge):
        for amplitude in (0.01, 0.7, 3.0):
            expected = source(amplitude*coefficients, **geometry)
            actual = source(amplitude*coefficients, **geometry, contiguous_arrays=True)
            relative_error = np.linalg.norm(actual-expected)/np.linalg.norm(expected)
            assert relative_error < 2e-12, relative_error
            print(f'Random multi-channel source, scale={amplitude}: '
                  f'relative L2 difference={relative_error:.3e}', flush=True)
        outputs = []
        for contiguous in (False, True):
            with patch.object(nonlinear_tower, 'nonlinear_source_on_grid',
                              partial(source, contiguous_arrays=contiguous)):
                outputs.append(coupled_rk4_step(towers, states, dt, **angular))
        for m in towers:
            for expected, actual in zip(outputs[0][m], outputs[1][m]):
                np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-14)
        print('Layout full-step agreement passed.', flush=True)
        elapsed = [[], []]
        for repeat in range(repeats):
            for index in ((0, 1) if repeat % 2 == 0 else (1, 0)):
                with patch.object(nonlinear_tower, 'nonlinear_source_on_grid',
                                  partial(source, contiguous_arrays=bool(index))):
                    start = perf_counter()
                    coupled_rk4_step(towers, states, dt, **angular)
                    elapsed[index].append(perf_counter()-start)
        for name, timings in zip(('Strided', 'Contiguous'), elapsed):
            print(f'{name}: median={np.median(timings):.6f}s; '
                  f'range={min(timings):.6f}..{max(timings):.6f}s', flush=True)
        print(f'Layout speedup: {np.median(elapsed[0])/np.median(elapsed[1]):.3f}x')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--julia', default='julia')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--check-tests', action='store_true')
    parser.add_argument('--profile-source', action='store_true',
                        help='time source operations separately using targeted line tracing')
    parser.add_argument('--compare-layouts', action='store_true',
                        help='compare strided and contiguous matrices with the Julia remainder')
    parser.add_argument('--profile', action='store_true',
                        help='profile three warmed Julia-backed RK4 steps after timing')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    numpy_remainder = nonlinear_source.sine_remainder
    startup = perf_counter()
    with JuliaRemainder(args.julia) as bridge:
        print(f'Julia startup/compile: {perf_counter()-startup:.3f}s (excluded from timings)', flush=True)
        # Exercise buffer reuse, scalar/strided/empty inputs and output ownership.
        values = np.linspace(-3, 3, 80).reshape(8, 10)[:, ::2]
        for field in (values, np.array(0.01), np.array([]), np.array([1e-100, -.1, .1])):
            saved = field.copy()
            result = bridge(field)
            np.testing.assert_allclose(result, numpy_remainder(field), rtol=5e-13, atol=0)
            snapshot = result.copy()
            bridge(np.ones(20))
            np.testing.assert_array_equal(result, snapshot)
            np.testing.assert_array_equal(field, saved)
        try:
            bridge(np.array([1j]))
        except ValueError:
            pass
        else:
            raise AssertionError('complex input must be rejected')
        if args.check_tests:
            with patch.object(nonlinear_source, 'sine_remainder', bridge):
                suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parents[2] / 'tests'))
                result = unittest.TextTestRunner().run(suite)
                if not result.wasSuccessful():
                    raise RuntimeError('Julia-backed regression tests failed')

        towers = {m: Tower(spin=.99, alpha=.42, m=m, l_max=5, x_left=-200,
                           x_right=511, grid_spacing_n=6, x_pml_center=-100,
                           pml_width=10) for m in (1, 3, 5)}
        reference = towers[1]
        modes = [(l, m) for m, tower in towers.items() for l in tower.l_values]
        cosine_theta, weights = leggauss(32)
        theta = np.arccos(cosine_theta)[:, None]
        azimuth = (2 * np.pi * np.arange(64) / 64)[None, :]
        angular = dict(harmonics=np.array([sph_harm_y(l, m, theta, azimuth) for l, m in modes]),
                       cosine_theta=cosine_theta, solid_angle_weights=weights[:, None]*2*np.pi/64)
        states = {m: tuple(np.zeros((tower.n_channels, tower.num_points), complex)
                           for _ in range(3)) for m, tower in towers.items()}
        profile = np.exp(-((reference.x-30)/10)**2)
        profile[:2] = profile[-2:] = 0
        peak = np.max(2*np.sqrt(3/(8*np.pi))*profile/reference.r)
        states[1][0][0] = profile/peak
        states[1][1][0] = -.4j*states[1][0][0]
        dt = .8*reference.dx
        if args.profile_source:
            profile_source_parts(towers, states, dt, angular, bridge, args.repeats)
            return
        if args.compare_layouts:
            compare_layouts(towers, states, dt, angular, bridge, args.repeats)
            return
        outputs = []
        for implementation in (numpy_remainder, bridge):
            with patch.object(nonlinear_source, 'sine_remainder', implementation):
                outputs.append(coupled_rk4_step(towers, states, dt, **angular))
        for m in towers:
            for expected, actual in zip(outputs[0][m], outputs[1][m]):
                np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-14)
        print('Full-step agreement passed; timing alternating NumPy/Julia steps.', flush=True)
        elapsed = [[], []]
        for repeat in range(args.repeats):
            for index in ((0, 1) if repeat % 2 == 0 else (1, 0)):
                with patch.object(nonlinear_source, 'sine_remainder', (numpy_remainder, bridge)[index]):
                    start = perf_counter()
                    coupled_rk4_step(towers, states, dt, **angular)
                    elapsed[index].append(perf_counter()-start)
        for name, timings in zip(('NumPy', 'Julia bridge'), elapsed):
            median = float(np.median(timings))
            print(f'{name}: median={median:.6f}s, range={min(timings):.6f}..{max(timings):.6f}s; '
                  f'1e5 M compute-only estimate={median*(1e5/dt)/3600:.2f}h')
        print(f'Full-step speedup: {np.median(elapsed[0])/np.median(elapsed[1]):.3f}x')
        if args.profile:
            profiler = cProfile.Profile()
            with patch.object(nonlinear_source, 'sine_remainder', bridge):
                for _ in range(3):
                    profiler.runcall(coupled_rk4_step, towers, states, dt, **angular)
            print('\nProfile totals for THREE warmed Julia-backed steps (not per step):')
            pstats.Stats(profiler).strip_dirs().sort_stats('cumulative').print_stats(22)


if __name__ == '__main__':
    main()
