"""Compare fresh stationary solves with the trajectory's interpolated inputs.

These checks distinguish time integration, interpolation, and radial/endpoint
accuracy. They do not test omitted self-gravity or higher cloud modes.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
from scipy.optimize import brentq

from adiabatic.primary_secondary.qcd_spindown import (
    INITIAL_ALPHA, INITIAL_CLOUD_FRACTIONS, StationaryTrack, evolve, physical_state,
)


RESULTS = Path(__file__).resolve().parents[1] / 'results'


def read_point(path):
    data = np.genfromtxt(path, names=True, delimiter=',')
    return {name: float(data[name]) for name in data.dtype.names}


def horizon_response(row, name, frequency, azimuthal_number):
    """Remove the vanishing horizon factor before comparing near saturation."""
    horizon_velocity = row['spin'] / (2 * (1 + np.sqrt(1 - row['spin']**2)))
    return row[name + '_horizon'] / (row[frequency] - azimuthal_number * horizon_velocity)


def main(stationary_directory, radial_directory, validation_directory):
    track = StationaryTrack(stationary_directory)
    trajectory = evolve(track)
    checks = []

    def compare(name, value, reference, tolerance):
        relative_error = abs(value / reference - 1)
        checks.append(dict(
            name=name, value=value, reference=reference,
            relative_error=relative_error, tolerance=tolerance,
            passed=bool(np.isfinite(relative_error) and relative_error < tolerance),
        ))

    def compare_rows(label, actual, reference, tolerance):
        for observable in ('omega1', 'omega2', 'm0_horizon', 'm3_outward'):
            compare(label + '_' + observable, actual[observable], reference[observable], tolerance)
        for name, frequency, azimuthal_number in (
            ('primary', 'omega1', 1), ('secondary', 'omega2', 2)
        ):
            # At the initial high spin, the much smaller l=2 horizon flux
            # changes by about 0.5% under radial refinement. Accept 1% for
            # that observable explicitly; do not relax the primary checks.
            flux_tolerance = 0.01 if name == 'secondary' else tolerance
            compare(
                label + '_' + name + '_horizon_response',
                horizon_response(actual, name, frequency, azimuthal_number),
                horizon_response(reference, name, frequency, azimuthal_number),
                flux_tolerance,
            )

    for index in (0, 11):
        actual = read_point(radial_directory / f'refine_point{index:03d}.csv')
        reference = read_point(stationary_directory / f'point{index:03d}.csv')
        compare_rows(f'radial_point{index:03d}', actual, reference, 0.003)

    for label, fraction in [('interior_one', .0275), ('interior_two', .0675)]:
        time = brentq(
            lambda time: INITIAL_CLOUD_FRACTIONS[0] * np.exp(trajectory.sol(time)[2]) - fraction,
            0, 1,
        )
        state = trajectory.sol(time)
        reference = track(state)
        reference['spin'] = physical_state(state)[1]
        actual = read_point(validation_directory / (label + '.csv'))
        compare_rows(label, actual, reference, 0.003)

    late_time = brentq(lambda time: physical_state(trajectory.sol(time))[1] - .52, 0, 1)
    state = trajectory.sol(late_time)
    reference = track(state)
    reference['spin'] = physical_state(state)[1]
    for label in ('late_refined', 'late_outer'):
        actual = read_point(validation_directory / (label + '.csv'))
        compare_rows(label, actual, reference, 0.003)

    # The independent Leaver solve uses depth 2000. A depth-500 solve at this
    # small alpha is not yet converged in the tiny imaginary frequency.
    initial = read_point(stationary_directory / 'point000.csv')
    compare('initial_Leaver_growth', -.5 * initial['primary_horizon'],
            1.8005587690408697e-10, 2e-5)
    compare('initial_Leaver_secondary_growth', -.5 * initial['secondary_horizon'],
            6.564018195700955e-16, 0.01)

    for index, row in enumerate(track.rows):
        for channel in ('m0', 'm3'):
            compare(f'point{index:03d}_{channel}_power', row[channel + '_power_ratio'], 1, 1e-3)

    with open(RESULTS / 'validation_checks.csv', 'w') as output:
        writer = csv.DictWriter(output, fieldnames=checks[0].keys())
        writer.writeheader()
        writer.writerows(checks)
    failures = [check for check in checks if not check['passed']]
    print(json.dumps(dict(total=len(checks), failures=failures), indent=2))
    if failures:
        raise RuntimeError('Stationary validation failed; inspect validation_checks.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stationary-directory', type=Path, default=RESULTS / 'stationary_v2')
    parser.add_argument('--radial-directory', type=Path, default=RESULTS / 'radial_checks')
    parser.add_argument('--validation-directory', type=Path, default=RESULTS / 'validation')
    arguments = parser.parse_args()
    main(arguments.stationary_directory, arguments.radial_directory, arguments.validation_directory)
