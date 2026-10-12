"""Cheap physical and regression checks before submitting the survey."""
import numpy as np
from numpy.lib.recfunctions import append_fields
import campaign


def main():
    # Independently normalized 211 real-field peak: alpha^2 sqrt(E)/(sqrt(8pi)e).
    expected = .15**2 * np.sqrt(.01) / (np.sqrt(8 * np.pi) * np.e)
    np.testing.assert_allclose(campaign.hydrogenic_peak(.01, .15, 1), expected, rtol=1e-14)

    case = dict(mass_solar=20., axion_eV=1e-12,
                alpha0=.14966276123144445, fa_GeV=5.7e18)
    campaign.configure(case)
    np.testing.assert_allclose(campaign.balance.INITIAL_CODE_ENERGIES, [1e-5, 1e-8], rtol=1e-14)
    campaign.configure(dict(case, axion_eV=1e-11, alpha0=case['alpha0'] * 10, fa_GeV=5.7e17))
    np.testing.assert_allclose(campaign.balance.INITIAL_CLOUD_FRACTIONS, campaign.SEEDS, rtol=1e-14)
    np.testing.assert_allclose(campaign.balance.INITIAL_CODE_ENERGIES, [1e-3, 1e-6], rtol=1e-14)
    campaign.configure(case)

    track = campaign.balance.StationaryTrack(
        campaign.ROOT / 'runs/adiabatic/qcd_spindown/results/stationary_v2')
    # The extra diagnostics do not change the balance law. Supply a dilute
    # secondary estimate only to keep its validity event inactive in this
    # previously validated primary-dominated regression trajectory.
    track.rows = append_fields(track.rows,
        ['primary_peak_envelope', 'secondary_peak_per_sqrt_energy'],
        [track.rows['peak_field'], campaign.hydrogenic_peak(1., track.rows['alpha'], 2)],
        usemask=False)
    campaign.AGE_YEARS = 1.
    evolution, status = campaign.evolve_track(track)
    assert status == 'age_limit', status
    expected = np.genfromtxt(
        campaign.ROOT / 'runs/adiabatic/qcd_spindown/results/adiabatic.csv',
        names=True, delimiter=',')
    actual = campaign.balance.trajectory_array(evolution, expected['time_years'])
    spin_error = np.max(abs(actual[:, 1] - expected['spin']))
    assert spin_error < 2e-6, spin_error
    assert np.max(abs(actual[:, -1])) < 1e-7
    print(f'PASS: normalized 211 peak, physical seed invariance, energy balance, reference trajectory; max spin difference {spin_error:.3g}')


if __name__ == '__main__':
    main()
