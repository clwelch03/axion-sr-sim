"""Restore the sampled primary field to GeV along the saved QCD trajectory.

Run with PYTHONPATH=src and the BHsuper Python. This is post-processing:
no stationary solves or time evolution are repeated. The measured field is
the real primary (including its conjugate and harmonics), sampled on the
stationary driver's grid; it is neither A0 nor the full two-cloud field.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from common import plotstyle


CAMPAIGN = Path(__file__).resolve().parents[1]
RESULTS = CAMPAIGN / 'results'
PLOTS = CAMPAIGN.parents[2] / 'plots'


def main():
    parameters_path = RESULTS / 'parameters.json'
    parameters = json.loads(parameters_path.read_text())
    decay_constant = parameters['fa_GeV']
    epsilon = parameters['epsilon']
    initial_alpha = parameters['alpha0']
    stationary_paths = sorted((RESULTS / 'stationary_v2').glob('point*.csv'))
    if len(stationary_paths) != 21:
        raise ValueError('Expected the complete 21-point reference table.')
    stationary = np.concatenate([
        np.atleast_1d(np.genfromtxt(path, names=True, delimiter=','))
        for path in stationary_paths
    ])
    trajectory_path = RESULTS / 'adiabatic.csv'
    trajectory = np.genfromtxt(trajectory_path, names=True, delimiter=',')

    # E_code = E_physical / (G fa^2 M_BH); fraction = E_physical / M_initial.
    # Remove the known sqrt(E_code) dependence before interpolating the slowly
    # varying profile shape. Directly splining peak_field across the very
    # unequal seed energies would misrepresent its square-root behavior.
    coordinates = epsilon * stationary['alpha'] / initial_alpha * stationary['energy']
    if np.any(np.diff(coordinates) <= 0):
        raise ValueError('Stationary energies must increase.')
    shape = stationary['peak_field'] / np.sqrt(stationary['energy'])
    shape_spline = CubicSpline(coordinates, shape, extrapolate=False)
    shape_pchip = PchipInterpolator(coordinates, shape, extrapolate=False)
    primary_fraction = trajectory['primary_fraction']
    if np.any(primary_fraction < coordinates[0]) or np.any(primary_fraction > coordinates[-1]):
        raise ValueError('Trajectory left the stationary table; refusing extrapolation.')
    mass_ratio = trajectory['mass_solar'] / parameters['mass_solar']
    code_energy = primary_fraction / (epsilon * mass_ratio)
    peak_dimensionless = shape_spline(primary_fraction) * np.sqrt(code_energy)
    peak_gev = decay_constant * peak_dimensionless
    cubic_over_mass = peak_dimensionless**2 / 6
    # Stable series for 1-sin(phi)/phi in this small-field run.
    sine_remainder_over_mass = cubic_over_mass * (
        1 - peak_dimensionless**2 / 20 + peak_dimensionless**4 / 840
    )
    frequency = CubicSpline(coordinates, stationary['omega1'], extrapolate=False)(primary_fraction)
    binding_scale = 1 - (frequency / trajectory['alpha'])**2
    if np.any(binding_scale <= 0):
        raise ValueError('Expected a bound primary frequency.')
    # A local scale estimate, not a computed frequency shift or an error bound
    # for perturbation theory: compare with (mu^2 - omega_R^2) Phi.
    cubic_over_binding = cubic_over_mass / binding_scale
    if not np.all(np.isfinite(peak_dimensionless)) or np.any(peak_dimensionless <= 0):
        raise ValueError('Invalid reconstructed field.')

    checks = []

    def check(name, value, reference, tolerance):
        relative_error = float(abs(value / reference - 1))
        checks.append(dict(name=name, relative_error=relative_error, tolerance=tolerance))
        if not np.isfinite(relative_error) or relative_error > tolerance:
            raise AssertionError(f'{name}: {relative_error} > {tolerance}')

    check('initial_field', peak_dimensionless[0], stationary['peak_field'][0], 1e-6)
    pchip_peak = shape_pchip(primary_fraction) * np.sqrt(code_energy)
    check('cubic_vs_pchip', 1 + np.max(abs(pchip_peak / peak_dimensionless - 1)), 1, 1e-4)
    validation_paths = sorted((RESULTS / 'validation').glob('*.csv'))
    if len(validation_paths) != 4:
        raise ValueError('Expected four independent validation solves.')
    for path in validation_paths:
        actual = np.genfromtxt(path, names=True, delimiter=',')
        coordinate = epsilon * actual['alpha'] / initial_alpha * actual['energy']
        prediction = float(shape_spline(coordinate) * np.sqrt(actual['energy']))
        check(path.stem + '_peak_field', prediction, actual['peak_field'], 1e-4)
    radial_paths = sorted((RESULTS / 'radial_checks').glob('*.csv'))
    if len(radial_paths) != 2:
        raise ValueError('Expected two radial refinement solves.')
    for path in radial_paths:
        actual = np.genfromtxt(path, names=True, delimiter=',')
        original = np.genfromtxt(
            RESULTS / 'stationary_v2' / path.name.removeprefix('refine_'),
            names=True, delimiter=',',
        )
        check(path.stem + '_peak_field', actual['peak_field'], original['peak_field'], 1e-4)

    columns = np.column_stack((
        trajectory['time_years'], trajectory['time_years'] * 365.25,
        peak_dimensionless, peak_gev, cubic_over_mass,
        sine_remainder_over_mass, cubic_over_binding,
    ))
    np.savetxt(RESULTS / 'field_amplitude.csv', columns, delimiter=',', comments='', header=(
        'time_years,time_days,primary_peak_Phi_over_fa,primary_peak_Phi_GeV,'
        'cubic_over_mass,sine_remainder_over_mass,cubic_over_binding_scale'
    ))

    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'stix', 'font.size': 11})
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    days = trajectory['time_years'] * 365.25
    axes[0].plot(days, peak_gev / 1e15, lw=2)
    axes[0].set_ylabel(r'Sampled primary peak $|\Phi|$ [$10^{15}$ GeV]')
    axes[0].text(.96, .12, rf'Final $|\Phi|/f_a = {peak_dimensionless[-1]:.2g}$',
                 ha='right', transform=axes[0].transAxes)
    axes[1].semilogy(days, cubic_over_mass, lw=2, label=r'Cubic / mass: $(\Phi/f_a)^2/6$')
    axes[1].semilogy(days, cubic_over_binding, lw=2, color=plotstyle.PALETTE[1],
                     label=r'Cubic / binding scale (estimate)')
    axes[1].set_ylabel('Dimensionless ratio')
    axes[1].legend(frameon=False, loc='lower right', fontsize=9)
    for axis in axes:
        axis.set_xlabel('Time from classical seed [days]')
        axis.set_xlim(0, 100)
        axis.grid(alpha=.2)
    figure.suptitle(
        r'QCD run: $M_0=20\,M_\odot$, $m_a=10^{-12}$ eV, $f_a=5.7\times10^{18}$ GeV'
        '\nPrimary field envelope; approximately constant from day 60 through one year',
        fontsize=12,
    )
    PLOTS.mkdir(exist_ok=True)
    for extension in ('png', 'pdf'):
        figure.savefig(PLOTS / f'qcd_field_amplitude_20Msun_1e-12eV.{extension}', dpi=180)
    plt.close(figure)

    summary = dict(
        fa_GeV=decay_constant,
        lambda_cosine=(parameters['axion_eV'] * 1e-9 / decay_constant)**2,
        initial_primary_peak_Phi_over_fa=float(peak_dimensionless[0]),
        final_primary_peak_Phi_over_fa=float(peak_dimensionless[-1]),
        initial_primary_peak_Phi_GeV=float(peak_gev[0]),
        final_primary_peak_Phi_GeV=float(peak_gev[-1]),
        field_growth_factor=float(peak_dimensionless[-1] / peak_dimensionless[0]),
        max_cubic_over_mass=float(np.max(cubic_over_mass)),
        max_cubic_over_binding_scale=float(np.max(cubic_over_binding)),
        stationary_table_max_peak=float(np.max(stationary['peak_field'])),
        final_secondary_to_primary_energy=float(
            trajectory['secondary_fraction'][-1] / primary_fraction[-1]),
        checks=checks,
        limitations=[
            'Sampled real primary field, including its higher harmonics; independent secondary and driven fields omitted.',
            'Interpolation checked against fresh solves, but angular sampling of the maximum is not refined.',
            'Binding-scale ratio is a local diagnostic, not a proof of convergence or a secular phase-error estimate.',
            'Inherits the preliminary scalar-only, cosine-potential, instantaneous-Kerr model of the saved trajectory.',
        ],
        input_sha256={
            str(path.relative_to(CAMPAIGN)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [parameters_path, trajectory_path, *stationary_paths, *validation_paths, *radial_paths]
        },
    )
    (RESULTS / 'field_amplitude_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    for name, value in summary.items():
        if name not in ('input_sha256', 'limitations', 'checks'):
            print(f'{name}: {value}')
    print(f'Passed {len(checks)} field checks; saved field_amplitude.csv, summary JSON, PNG and PDF.')


if __name__ == '__main__':
    main()
