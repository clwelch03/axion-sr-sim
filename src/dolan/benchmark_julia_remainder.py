"""Compare NumPy and Julia on identical float64 inputs; no Python-Julia bridge.

Run: PYTHONPATH=src python -m dolan.benchmark_julia_remainder --julia /path/to/julia
Requires Julia Base only. Temporary binary files preserve inputs exactly.
"""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile
from time import perf_counter

import mpmath as mp
import numpy as np

from dolan.nonlinear_source import sine_remainder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--julia', default='julia')
    args = parser.parse_args()
    random = np.random.default_rng(42)
    count = 2048 * 256
    # A real m=1 Gaussian radial slice, normalized to sampled peak Phi=0.7.
    u, _ = np.polynomial.legendre.leggauss(32)
    azimuth = 2 * np.pi * np.arange(64) / 64
    angular = (np.sqrt(1 - u[:, None]**2) * np.cos(azimuth)).ravel()
    radii = np.linspace(3, 100, 256)
    cloud = angular[:, None] * (np.exp(-((radii - 30) / 10)**2) / radii)
    cloud *= 0.7 / np.max(np.abs(cloud))
    cases = dict(small=random.uniform(-.099, .099, count),
                 mixed=random.uniform(-3, 3, count),
                 tiny=random.choice([-1., 1.], count) * 10**random.uniform(-160, -10, count),
                 cloud=cloud.ravel(),
                 edges=np.array([-3., -.1001, -.1, -.0999, -1e-8, -0., 0.,
                                 1e-160, 1e-100, 1e-12, .0999, .1, .1001, 3.]))
    timings = {}
    def time_numpy(name, field):
        sine_remainder(field)
        for _ in range(11):
            start = perf_counter()
            sine_remainder(field)
            timings.setdefault(name, []).append(perf_counter() - start)

    with tempfile.TemporaryDirectory(prefix='axion-julia-') as directory:
        for name, field in cases.items():
            field.tofile(Path(directory) / (name + '.input'))
            time_numpy(name, field)
        result = subprocess.run(
            [args.julia, '--startup-file=no', '--history-file=no', '--threads=1',
             str(Path(__file__).with_suffix('.jl').resolve()), directory],
            env=dict(os.environ, JULIA_NUM_THREADS='1'), check=True, capture_output=True, text=True)
        records = {line.split()[0]: list(map(float, line.split()[1:]))
                   for line in result.stdout.splitlines()}
        print('case     NumPy ms  Julia alloc ms  Julia reuse ms  max abs difference  Julia bytes (alloc/reuse)')
        for name, field in cases.items():
            time_numpy(name, field)  # Bracket Julia execution with NumPy timings.
            actual = np.fromfile(Path(directory) / (name + '.output'), dtype=np.float64)
            expected = sine_remainder(field)
            np.testing.assert_allclose(actual, expected, rtol=5e-13,
                                       atol=8 * np.nextafter(0., 1.))
            if name == 'edges':
                with mp.workdps(550):
                    accurate = np.array([float(mp.mpf(float(x)) - mp.sin(mp.mpf(float(x))))
                                         for x in field])
                np.testing.assert_allclose(actual, accurate, rtol=5e-13,
                                           atol=8 * np.nextafter(0., 1.))
            allocating, inplace, allocated_bytes, inplace_bytes = records[name]
            print(f'{name:8s} {1e3*np.median(timings[name]):9.3f} {1e3*allocating:15.3f}'
                  f' {1e3*inplace:15.3f} {np.max(np.abs(actual-expected)):19.3e}'
                  f' {int(allocated_bytes)}/{int(inplace_bytes)}')
        print('Accuracy checks passed. Kernel timings exclude compilation, I/O, and interoperability overhead.')


if __name__ == '__main__':
    main()
