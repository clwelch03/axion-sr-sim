"""Integrate the saved stationary track, audit it, and draw the comparison.

Run from the repository root with PYTHONPATH=src and the BHsuper Python.
The stationary CSV files are produced by submit.sh, not by this script.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from common import plotstyle
from adiabatic.primary_secondary.qcd_spindown import (
    CSV_HEADER, EPSILON, INITIAL_CLOUD_FRACTIONS, StationaryTrack,
    evolve, independent_coefficients, physical_state, save_metadata,
    trajectory_array,
)


CAMPAIGN = Path(__file__).resolve().parents[1]
RESULTS = CAMPAIGN / 'results'
PLOTS = CAMPAIGN.parents[2] / 'plots'


def time_at_spin(evolution, target_spin):
    """First threshold crossing during the monotonically decreasing episode."""
    return brentq(
        lambda time: physical_state(evolution.sol(time))[1] - target_spin,
        0, evolution.t[-1], xtol=1e-12,
    )


def main(stationary_directory):
    track = StationaryTrack(stationary_directory)
    if len(track.rows) != 21:
        raise RuntimeError('The full 21-point stationary sequence is required.')

    adiabatic = evolve(track)
    refined_time = evolve(track, rtol=2e-10, max_step=.0005)

    def nonrelativistic_coefficients(state):
        _, spin, alpha, _ = physical_state(state)
        return independent_coefficients(alpha, spin)

    nonrelativistic = evolve(nonrelativistic_coefficients)
    times = np.linspace(0, 1, 4001)
    adiabatic_values = trajectory_array(adiabatic, times)
    refined_values = trajectory_array(refined_time, times)
    nonrelativistic_values = trajectory_array(nonrelativistic, times)

    for name, values in (
        ('adiabatic', adiabatic_values),
        ('adiabatic_refined_time', refined_values),
        ('nonrelativistic', nonrelativistic_values),
    ):
        np.savetxt(
            RESULTS / (name + '.csv'), values, delimiter=',',
            header=CSV_HEADER, comments='',
        )

    comparison = np.column_stack([
        times, adiabatic_values[:, 1], nonrelativistic_values[:, 1],
        adiabatic_values[:, 2], nonrelativistic_values[:, 2],
        adiabatic_values[:, 4:6], nonrelativistic_values[:, 4:6],
    ])
    np.savetxt(
        RESULTS / 'spindown_comparison.csv', comparison, delimiter=',', comments='',
        header=(
            'time_years,spin_adiabatic,spin_nonrelativistic,'
            'mass_solar_adiabatic,mass_solar_nonrelativistic,'
            'primary_fraction_adiabatic,secondary_fraction_adiabatic,'
            'primary_fraction_nonrelativistic,secondary_fraction_nonrelativistic'
        ),
    )

    # Validate on the actual trajectory, including points between table nodes.
    # A radial refinement at the initial point separately checks discretization.
    validation_targets = []
    for name, fraction in [('interior_one', .0275), ('interior_two', .0675)]:
        crossing = brentq(
            lambda time: INITIAL_CLOUD_FRACTIONS[0] *
            np.exp(adiabatic.sol(time)[2]) - fraction,
            0, 1,
        )
        validation_targets.append((name, crossing, .25, 2000.))
    # Check near saturation while the horizon factor is still resolved.
    # At the final numerical fixed point both flux and threshold factor can
    # be at round-off, making their ratio a poor discretization diagnostic.
    late_time = time_at_spin(adiabatic, .52)
    validation_targets.extend([
        ('late_refined', late_time, .125, 2000.),
        ('late_outer', late_time, .25, 2400.),
    ])

    with open(RESULTS / 'validation_targets.csv', 'w') as output:
        output.write('label,alpha,spin,energy,step,outer\n')
        for name, time, step, outer in validation_targets:
            _, spin, alpha, energies = physical_state(adiabatic.sol(time))
            output.write(
                f'{name},{alpha:.17g},{spin:.17g},{energies[0]:.17g},'
                f'{step},{outer}\n'
            )

    alpha_offsets = []
    spin_offsets = []
    for state in adiabatic.sol(times).T:
        _, spin, alpha, _ = physical_state(state)
        coordinate = INITIAL_CLOUD_FRACTIONS[0] * np.exp(state[2])
        alpha_offsets.append(alpha - track.splines['alpha'](coordinate))
        spin_offsets.append(spin - track.splines['spin'](coordinate))

    diagnostics = dict(
        max_energy_balance_error=float(np.max(np.abs(adiabatic_values[:, -1]))),
        max_spin_change_under_time_refinement=float(
            np.max(np.abs(adiabatic_values[:, 1] - refined_values[:, 1]))
        ),
        max_alpha_offset_from_stationary_track=float(np.max(np.abs(alpha_offsets))),
        max_spin_offset_from_stationary_track=float(np.max(np.abs(spin_offsets))),
        max_primary_field_over_fa=float(np.max(track.rows['peak_field'])),
        adiabatic_time_to_spin075_years=time_at_spin(adiabatic, .75),
        nr_time_to_spin075_years=time_at_spin(nonrelativistic, .75),
        adiabatic_final_spin=float(adiabatic_values[-1, 1]),
        nr_final_spin=float(nonrelativistic_values[-1, 1]),
        adiabatic_final_mass_solar=float(adiabatic_values[-1, 2]),
        adiabatic_final_secondary_fraction=float(adiabatic_values[-1, 5]),
    )
    assert diagnostics['max_energy_balance_error'] < 2e-8
    assert diagnostics['max_spin_change_under_time_refinement'] < 2e-6
    # Keep the scientific inputs alongside code/data hashes, rather than
    # relying on filenames alone to identify which stationary table was used.
    source_snapshot = RESULTS / 'stationary_driver_v2_snapshot.jl'
    provenance = dict(
        repository_revision=subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=CAMPAIGN, text=True,
        ).strip(),
        stationary_source_sha256=hashlib.sha256(source_snapshot.read_bytes()).hexdigest(),
        stationary_csv_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(stationary_directory.glob('point*.csv'))
        },
    )
    save_metadata(RESULTS / 'parameters.json', diagnostics=diagnostics, provenance=provenance)
    print(json.dumps(diagnostics, indent=2))

    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    figure, (spin_axis, cloud_axis) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={'height_ratios': [2, 1]}, constrained_layout=True,
    )
    spin_axis.plot(times, adiabatic_values[:, 1], lw=2.4, label='Adiabatic: nonlinear 211 + linear 322')
    spin_axis.plot(times, nonrelativistic_values[:, 1], '--', lw=2.1, label='Linear nonrelativistic: 211 + 322')
    spin_axis.set_ylabel(r'BH spin $\chi=J/M^2$')
    spin_axis.set_ylim(.47, 1.01)
    spin_axis.legend(loc='upper right', fontsize=10)
    spin_axis.grid(alpha=.2)
    spin_axis.set_title(
        r'$M_0=20\,M_\odot$, $\chi_0=0.99$, $m_a=10^{-12}\,\mathrm{eV}$'
        '\nQCD coupling; leading scalar-only BH/cloud balance',
        fontsize=13,
    )
    cloud_axis.semilogy(times, adiabatic_values[:, 4], lw=2, label=r'$E_{211}/M_0$ (adiabatic)')
    cloud_axis.semilogy(times, adiabatic_values[:, 5], lw=2, color=plotstyle.PALETTE[2],
                        label=r'$E_{322}/M_0$ (adiabatic)')
    cloud_axis.set_ylabel('Cloud energy fraction')
    cloud_axis.set_xlabel('Time from the stated classical seeds [years]')
    cloud_axis.set_xlim(0, .4)
    cloud_axis.set_ylim(5e-10, .2)
    cloud_axis.grid(alpha=.2)
    cloud_axis.legend(loc='center right', fontsize=10)
    figure.supxlabel(
        r'$E_{211,0}/(Gf_a^2 M_0)=10^{-5}$; $E_{322,0}/(Gf_a^2 M_0)=10^{-8}$'
        '\nPreliminary: no GW losses, self-gravity, or higher independent clouds.',
        fontsize=9,
    )
    PLOTS.mkdir(exist_ok=True)
    for extension in ('png', 'pdf'):
        figure.savefig(PLOTS / ('qcd_spindown_20Msun_1e-12eV.' + extension), dpi=180)
    plt.close(figure)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--stationary-directory', type=Path, default=RESULTS / 'stationary_v2',
        help='Directory containing the 21 accepted stationary point CSVs.',
    )
    main(parser.parse_args().stationary_directory)
