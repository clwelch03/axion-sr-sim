#!/usr/bin/env python
"""Independent-channel linear driver for the multilevel scaffold.

This is deliberately the smallest useful multimode evolution layer.  Each
active ``ModeSpec`` gets its own linear radial operator and shared tortoise-grid
PML; interlevel couplings are not included yet.  With the active set reduced to
``[211]``, the driver is required to reproduce the existing l=m=1 tower.
"""

import os
import sys

import numpy as np


_HERE = os.path.dirname(os.path.abspath(__file__))
_RADIAL1D = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _RADIAL1D)

from radial1d import build_grid, coeffs, d_dx, d2_dx2, horizons, tortoise_x
from mode_catalog import ModeSpec, leaver_frequency
from radial_modes import build_radial_mode


def run(active_modes, alpha=0.30, dimensionless_spin=0.90,
        xmin=-40.0, xmax=120.0, grid_points=600, final_time=60.0,
        cfl=0.4, eps_match=0.02, pml_peak_damping=1.0,
        probe_radius=None, recorded_points=120):
    """Evolve independent linear channels for an explicit active mode set.

    The independent-channel structure is intentional: it provides the
    multimode data layout and reduction test before any quartic transfer terms
    are introduced.
    """
    active_modes = tuple(active_modes)
    if not active_modes:
        raise ValueError("active_modes must contain at least one ModeSpec")
    if not all(isinstance(mode_spec, ModeSpec) for mode_spec in active_modes):
        raise TypeError("active_modes must contain only ModeSpec instances")

    x, radius = build_grid(
        dimensionless_spin, xmin, xmax, grid_points
    )
    grid_spacing = x[1] - x[0]
    horizon_radius, _ = horizons(dimensionless_spin)
    horizon_angular_frequency = (
        dimensionless_spin / (2.0 * horizon_radius)
    )
    channel_count = len(active_modes)
    complex_frequencies = np.array([
        leaver_frequency(mode_spec, alpha, dimensionless_spin)
        for mode_spec in active_modes
    ])

    denominator = np.empty((channel_count, grid_points))
    second_derivative_coefficient = np.empty((channel_count, grid_points))
    first_derivative_coefficient = np.empty((channel_count, grid_points), complex)
    time_derivative_coefficient = np.empty((channel_count, grid_points), complex)
    potential = np.empty((channel_count, grid_points), complex)
    for channel_index, mode_spec in enumerate(active_modes):
        diagonal_cosine_coupling = (
            1.0 / 3.0
            + (2.0 / 3.0)
            * (mode_spec.ell * (mode_spec.ell + 1.0)
               - 3.0 * mode_spec.azimuthal_number**2)
            / ((2.0 * mode_spec.ell - 1.0)
               * (2.0 * mode_spec.ell + 3.0))
        )
        coefficients = coeffs(
            radius, dimensionless_spin, alpha,
            l=mode_spec.ell, m=mode_spec.azimuthal_number,
            c_ll=diagonal_cosine_coupling
        )
        (denominator[channel_index],
         second_derivative_coefficient[channel_index],
         first_derivative_coefficient[channel_index],
         time_derivative_coefficient[channel_index],
         potential[channel_index], _, _) = coefficients
    inverse_denominator = 1.0 / denominator

    x_match = tortoise_x(horizon_radius + eps_match, dimensionless_spin)
    pml_mask = x < x_match
    attenuation = np.zeros_like(x)
    attenuation[pml_mask] = pml_peak_damping * (
        (x_match - x[pml_mask]) / (x_match - xmin)
    ) ** 2
    outer_sponge_start = xmax - 40.0
    outer_sponge_mask = x > outer_sponge_start
    outer_damping = np.zeros_like(x)
    outer_damping[outer_sponge_mask] = 0.8 * (
        (x[outer_sponge_mask] - outer_sponge_start) / 40.0
    ) ** 2
    physical_mask = ~pml_mask & ~outer_sponge_mask
    imaginary_horizon_frequencies = np.array([
        1j * mode_spec.azimuthal_number * horizon_angular_frequency
        for mode_spec in active_modes
    ])

    if probe_radius is None:
        probe_radius = min(2.0 / alpha**2, 0.6 * radius[-1])
    probe_index = int(np.argmin(np.abs(radius - probe_radius)))

    time_step = cfl * grid_spacing
    step_count = int(final_time / time_step)
    recorded_times = np.linspace(0.0, final_time, recorded_points)
    recording_indices = (recorded_times / time_step).astype(int)
    next_recording = 0

    reduced_field = np.zeros((channel_count, grid_points), complex)
    normal_momentum = np.zeros_like(reduced_field)
    pml_auxiliary = np.zeros_like(reduced_field)
    for channel_index, mode_spec in enumerate(active_modes):
        radial_mode = build_radial_mode(
            mode_spec, alpha, dimensionless_spin, radius
        )
        reduced_field[channel_index] = radial_mode.reduced_radial_profile
        normal_momentum[channel_index] = radial_mode.reduced_momentum_profile

    def right_hand_side(field, momentum, auxiliary):
        field_time_derivative = np.empty_like(field)
        momentum_time_derivative = np.empty_like(momentum)
        auxiliary_time_derivative = np.empty_like(auxiliary)
        for channel_index in range(channel_count):
            field_x = d_dx(field[channel_index], grid_spacing)
            field_xx = d2_dx2(field[channel_index], grid_spacing)
            auxiliary_x = d_dx(auxiliary[channel_index], grid_spacing)
            interior_momentum_rhs = inverse_denominator[channel_index] * (
                second_derivative_coefficient[channel_index] * field_xx
                + first_derivative_coefficient[channel_index] * field_x
                + time_derivative_coefficient[channel_index]
                * momentum[channel_index]
                - potential[channel_index] * field[channel_index]
            )
            pml_momentum_rhs = (
                field_xx
                + 2.0 * imaginary_horizon_frequencies[channel_index]
                * (field_x - momentum[channel_index])
                - auxiliary_x
                - imaginary_horizon_frequencies[channel_index]
                * auxiliary[channel_index]
                - attenuation * (
                    momentum[channel_index]
                    + imaginary_horizon_frequencies[channel_index]
                    * field[channel_index]
                )
            )
            momentum_time_derivative[channel_index] = (
                np.where(pml_mask, pml_momentum_rhs, interior_momentum_rhs)
                - outer_damping * momentum[channel_index]
            )
            field_time_derivative[channel_index] = (
                momentum[channel_index] - outer_damping * field[channel_index]
            )
            auxiliary_time_derivative[channel_index] = np.where(
                pml_mask,
                attenuation * (
                    field_x - auxiliary[channel_index]
                    + imaginary_horizon_frequencies[channel_index]
                    * field[channel_index]
                )
                - imaginary_horizon_frequencies[channel_index]
                * auxiliary[channel_index],
                0.0
            )
        return (field_time_derivative, momentum_time_derivative,
                auxiliary_time_derivative)

    recorded_field = []
    recorded_cloud_norm = []
    recorded_time_values = []
    for step_index in range(step_count):
        if (next_recording < recorded_points
                and step_index >= recording_indices[next_recording]):
            recorded_time_values.append(step_index * time_step)
            recorded_field.append(reduced_field[:, probe_index].copy())
            recorded_cloud_norm.append(np.sqrt(
                np.sum(np.abs(reduced_field[:, physical_mask])**2, axis=1)
                * grid_spacing
            ))
            next_recording += 1

        first_stage = right_hand_side(
            reduced_field, normal_momentum, pml_auxiliary
        )
        second_stage = right_hand_side(
            reduced_field + 0.5 * time_step * first_stage[0],
            normal_momentum + 0.5 * time_step * first_stage[1],
            pml_auxiliary + 0.5 * time_step * first_stage[2]
        )
        third_stage = right_hand_side(
            reduced_field + 0.5 * time_step * second_stage[0],
            normal_momentum + 0.5 * time_step * second_stage[1],
            pml_auxiliary + 0.5 * time_step * second_stage[2]
        )
        fourth_stage = right_hand_side(
            reduced_field + time_step * third_stage[0],
            normal_momentum + time_step * third_stage[1],
            pml_auxiliary + time_step * third_stage[2]
        )
        reduced_field += (time_step / 6.0) * (
            first_stage[0] + 2.0 * second_stage[0]
            + 2.0 * third_stage[0] + fourth_stage[0]
        )
        normal_momentum += (time_step / 6.0) * (
            first_stage[1] + 2.0 * second_stage[1]
            + 2.0 * third_stage[1] + fourth_stage[1]
        )
        pml_auxiliary += (time_step / 6.0) * (
            first_stage[2] + 2.0 * second_stage[2]
            + 2.0 * third_stage[2] + fourth_stage[2]
        )
        reduced_field[:, 0] = 0.0
        normal_momentum[:, 0] = 0.0
        reduced_field[:, -1] = 0.0
        normal_momentum[:, -1] = 0.0
        if not np.all(np.isfinite(reduced_field)):
            raise FloatingPointError(f"non-finite field at t={step_index * time_step}")

    return {
        "active_modes": active_modes,
        "complex_frequencies": complex_frequencies,
        "t": np.asarray(recorded_time_values),
        "probe_field": np.asarray(recorded_field).T,
        "cloud_norm": np.asarray(recorded_cloud_norm).T,
        "probe_radius": radius[probe_index],
        "grid_spacing": grid_spacing,
        "time_step": time_step,
        "x": x,
        "radius": radius,
    }


if __name__ == "__main__":
    from mode_catalog import MODE_211

    output = run((MODE_211,))
    print("active modes:", [mode.label for mode in output["active_modes"]])
    print("recorded field shape:", output["probe_field"].shape)
