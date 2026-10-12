"""Leading instantaneous-background BH/two-cloud evolution in fixed M0 units.

The Julia solver supplies nonlinear 211, averaged-background 322 and driven
m=0,3 fluxes. This module only interpolates those inputs and integrates balance.
The nonrelativistic comparison uses two independent linear hydrogenic clouds.
"""
from pathlib import Path
import json

import numpy as np
from scipy.constants import G, c, hbar, electron_volt
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from common.hydrogenic import growth_rate

SOLAR_MASS_KG = 1.988409870698051e30  # IAU nominal GM_sun divided by G
YEAR_SECONDS = 365.25 * 86400
MASS_SOLAR = 20.0
AXION_EV = 1e-12
INITIAL_SPIN = 0.99
DECAY_CONSTANT_GEV = 5.7e18
REFERENCE_TIME = G * MASS_SOLAR * SOLAR_MASS_KG / c**3
INITIAL_ALPHA = REFERENCE_TIME * AXION_EV * electron_volt / hbar
PLANCK_ENERGY_GEV = np.sqrt(hbar * c**5 / G) / (1e9 * electron_volt)
EPSILON = (DECAY_CONSTANT_GEV / PLANCK_ENERGY_GEV)**2
INITIAL_CODE_ENERGIES = np.array([1e-5, 1e-8])
INITIAL_CLOUD_FRACTIONS = EPSILON * INITIAL_CODE_ENERGIES


def physical_state(state):
    """Convert the evolved state to current BH units.

    state = [1-b, chi0-j, log(q1/q1_initial), log(q2/q2_initial), E_out/M0].
    Here b=M_BH/M0 and j=J_BH/M0^2. Log energies preserve positivity and
    resolve both the growing primary and the much smaller secondary cloud.
    """
    mass_ratio = 1 - state[0]
    angular_momentum = INITIAL_SPIN - state[1]
    cloud_fractions = INITIAL_CLOUD_FRACTIONS * np.exp(state[2:4])
    return (
        mass_ratio,
        angular_momentum / mass_ratio**2,
        INITIAL_ALPHA * mass_ratio,
        cloud_fractions / (EPSILON * mass_ratio),
    )


def independent_coefficients(alpha, spin):
    """Signed small-alpha rates; omega_R=mu matches the rate's threshold."""
    return dict(
        omega1=alpha,
        omega2=alpha,
        secondary_slope=0.0,
        primary_horizon=-2 * float(growth_rate(alpha, spin)),
        primary_outward=0.0,
        secondary_horizon=-2 * float(growth_rate(alpha, spin, l=2, m=2)),
        secondary_outward=0.0,
        m0_horizon=0.0,
        m0_outward=0.0,
        m3_horizon=0.0,
        m3_outward=0.0,
    )


def coupled_rhs(state, coefficients):
    """Derivative per year, identical to coupled_balance.jl at leading order.

    Cloud energies q_i use initial BH units. Flux coefficients use current BH
    units. The powers below restore each channel's actual normalization.
    Channel-by-channel projection avoids subtracting the large primary terms
    to obtain the extremely small secondary rate.
    """
    mass_ratio, _, _, energies = physical_state(state)
    primary_energy, secondary_energy = energies
    omega1, omega2 = coefficients['omega1'], coefficients['omega2']
    secondary_response = (
        2 * secondary_energy * coefficients['secondary_slope'] / omega2**2
    )
    primary_angular_derivative = 1 / omega1 - secondary_response
    secondary_angular_derivative = 2 / omega2
    determinant = secondary_angular_derivative - primary_angular_derivative

    primary_rate = secondary_rate = 0.0
    horizon_energy = horizon_torque = outward_energy = 0.0
    channels = [
        ('primary', primary_energy, 1 / omega1),
        ('secondary', secondary_energy, 2 / omega2),
        ('m0', primary_energy**2 * secondary_energy, 0.0),
        ('m3', primary_energy * secondary_energy**2, 3 / (2 * omega2 - omega1)),
    ]
    for name, normalization, torque_ratio in channels:
        horizon = coefficients[name + '_horizon'] * normalization
        outward = coefficients[name + '_outward'] * normalization
        total = horizon + outward
        primary_rate += total * (torque_ratio - secondary_angular_derivative) / determinant

        # Preserve the small response term when the channel is the primary.
        angular_difference = (
            -secondary_response if name == 'primary'
            else primary_angular_derivative - torque_ratio
        )
        secondary_rate += total * angular_difference / determinant
        horizon_energy += horizon
        horizon_torque += torque_ratio * horizon
        outward_energy += outward

    derivative = np.array([
        -EPSILON * horizon_energy,
        -EPSILON * mass_ratio * horizon_torque,
        primary_rate / (mass_ratio * primary_energy),
        secondary_rate / (mass_ratio * secondary_energy),
        EPSILON * outward_energy,
    ])
    return derivative * YEAR_SECONDS / REFERENCE_TIME


def evolve(coefficient_function, final_years=1.0, rtol=2e-8, max_step=0.002):
    """Adaptive Runge–Kutta evolution with an accumulated outgoing-energy audit."""
    def rhs(time, state):
        mass_ratio, spin, _, _ = physical_state(state)
        if not (mass_ratio > 0 and 0 < spin < 1):
            raise ValueError('left positive-mass, subextremal Kerr domain')
        return coupled_rhs(state, coefficient_function(state))

    result = solve_ivp(
        rhs, (0, final_years), np.zeros(5), rtol=rtol,
        atol=[1e-12, 1e-11, 1e-9, 1e-9, 1e-16],
        max_step=max_step, dense_output=True,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result


def write_targets(path):
    """Pilot path uses binding-corrected frequencies only to choose solve sites.

    It is not the adiabatic result. Subsequent evolution must verify that its
    alpha/spin track remains close, and fresh stationary checks test that error.
    """
    # On the pilot ignore the tiny secondary when relating extracted mass to J.
    def angular_derivative(extracted, angular_momentum):
        mass_ratio = 1 - extracted
        alpha = INITIAL_ALPHA * mass_ratio
        return [-mass_ratio / (alpha * (1 - alpha**2 / 8))]

    pilot = solve_ivp(
        angular_derivative, (0, 0.09), [INITIAL_SPIN],
        rtol=1e-12, atol=1e-13, dense_output=True,
    )
    losses = np.r_[0.0, 1e-5, 1e-4, 0.001, np.linspace(0.005, 0.085, 17)]
    with open(path, 'w') as output:
        output.write('label,alpha,spin,energy,step,outer\n')
        for index, extracted in enumerate(losses):
            mass_ratio = 1 - extracted
            spin = pilot.sol(extracted)[0] / mass_ratio**2
            energy = (INITIAL_CLOUD_FRACTIONS[0] + extracted) / (EPSILON * mass_ratio)
            output.write(
                f'point{index:03d},{INITIAL_ALPHA * mass_ratio:.17g},{spin:.17g},'
                f'{energy:.17g},0.25,2000\n'
            )


class StationaryTrack:
    """Cubic interpolation of dimensionally factored stationary fluxes.

    Its one-dimensional coordinate is x1=E211_physical/M0, including the seed.
    This is a trajectory-local surrogate, not a general alpha/spin/energy grid.
    """
    fields = (
        'omega1', 'omega2', 'secondary_slope',
        'primary_horizon', 'primary_outward',
        'secondary_horizon', 'secondary_outward',
        'm0_horizon', 'm0_outward', 'm3_horizon', 'm3_outward',
    )

    def __init__(self, directory):
        paths = sorted(Path(directory).glob('point*.csv'))
        self.rows = np.concatenate([
            np.atleast_1d(np.genfromtxt(path, names=True, delimiter=',')) for path in paths
        ])
        mass_ratios = self.rows['alpha'] / INITIAL_ALPHA
        self.coordinates = EPSILON * mass_ratios * self.rows['energy']
        self.splines = {
            name: CubicSpline(self.coordinates, self.rows[name], extrapolate=False)
            for name in self.fields + ('alpha', 'spin')
        }

        # Factor out the horizon superradiance threshold before interpolating.
        # The evolving spin need not equal the pilot spin exactly. Restoring
        # its CURRENT horizon factor prevents a small pilot error from moving
        # the saturation point, or giving the flux the wrong sign there.
        horizon_velocity = self.rows['spin'] / (
            2 * (1 + np.sqrt(1 - self.rows['spin']**2))
        )
        for name, frequency, azimuthal_number in (
            ('primary', 'omega1', 1), ('secondary', 'omega2', 2)
        ):
            threshold_factor = self.rows[frequency] - azimuthal_number * horizon_velocity
            response = self.rows[name + '_horizon'] / threshold_factor
            self.splines[name + '_horizon_response'] = CubicSpline(
                self.coordinates, response, extrapolate=False
            )

    def __call__(self, state):
        coordinate = INITIAL_CLOUD_FRACTIONS[0] * np.exp(state[2])
        if not self.coordinates[0] * (1 - 1e-6) <= coordinate <= self.coordinates[-1]:
            raise ValueError(f'primary cloud outside stationary table: {coordinate}')
        coordinate = max(coordinate, self.coordinates[0])
        result = {name: float(self.splines[name](coordinate)) for name in self.fields}

        _, spin, alpha, _ = physical_state(state)
        pilot_spin = float(self.splines['spin'](coordinate))
        pilot_alpha = float(self.splines['alpha'](coordinate))
        if abs(spin - pilot_spin) > 0.002 or abs(alpha / pilot_alpha - 1) > 1e-4:
            raise ValueError('trajectory left the stationary track; make fresh solves')

        horizon_velocity = spin / (2 * (1 + np.sqrt(1 - spin**2)))
        for name, frequency, azimuthal_number in (
            ('primary', 'omega1', 1), ('secondary', 'omega2', 2)
        ):
            response = float(self.splines[name + '_horizon_response'](coordinate))
            result[name + '_horizon'] = response * (
                result[frequency] - azimuthal_number * horizon_velocity
            )
        return result


def trajectory_array(evolution, times):
    """Export physical observables and the integrated energy-balance residual."""
    states = evolution.sol(times)
    mass_ratio = 1 - states[0]
    spin = (INITIAL_SPIN - states[1]) / mass_ratio**2
    fractions = INITIAL_CLOUD_FRACTIONS[:, None] * np.exp(states[2:4])
    energy_audit = (
        -states[0]
        + np.sum(fractions - INITIAL_CLOUD_FRACTIONS[:, None], axis=0)
        + states[4]
    )
    return np.column_stack([
        times, spin, MASS_SOLAR * mass_ratio, INITIAL_ALPHA * mass_ratio,
        fractions.T, states[4], energy_audit,
    ])


CSV_HEADER = (
    'time_years,spin,mass_solar,alpha,primary_fraction,secondary_fraction,'
    'emitted_fraction,energy_audit'
)


def save_metadata(path, **extra):
    metadata = dict(
        mass_solar=MASS_SOLAR, axion_eV=AXION_EV, spin0=INITIAL_SPIN,
        fa_GeV=DECAY_CONSTANT_GEV, alpha0=INITIAL_ALPHA, epsilon=EPSILON,
        seed_code_energies=INITIAL_CODE_ENERGIES.tolist(),
        seed_mass_fractions=INITIAL_CLOUD_FRACTIONS.tolist(),
        **extra,
    )
    Path(path).write_text(json.dumps(metadata, indent=2) + '\n')
