"""Dolan's frequency filter for isolating one time-domain mode.

Reference: Dolan, Phys. Rev. D 87, 124026 (2013), arXiv:1212.1477,
Sec. IV A 2, Eq. (49).
"""

import numpy as np


def fourier_amplitudes(times, signal):
    """Evaluate Dolan Eq. (46) on a uniformly sampled complex time series.

    The returned angular frequencies are signed and ordered from negative to
    positive.  The transform convention is exp(-i omega t), exactly as in
    Eq. (46).  Consequently, a solver mode with time dependence
    exp(-i omega_R t) produces a peak at Fourier frequency -omega_R.
    """
    times = np.asarray(times, dtype=float)
    signal = np.asarray(signal)
    if times.ndim != 1 or signal.ndim != 1:
        raise ValueError("times and signal must both be one-dimensional")
    if times.size != signal.size:
        raise ValueError("times and signal must have the same length")
    if times.size < 2:
        raise ValueError("at least two time samples are required")

    time_steps = np.diff(times)
    sample_spacing = time_steps[0]
    if sample_spacing <= 0.0 or not np.allclose(
            time_steps, sample_spacing, rtol=1.0e-10,
            atol=1.0e-12 * sample_spacing):
        raise ValueError("times must be strictly increasing and uniformly sampled")

    angular_frequencies = 2.0 * np.pi * np.fft.fftfreq(
        times.size, d=sample_spacing)
    amplitudes = np.fft.fft(signal) / times.size
    return np.fft.fftshift(angular_frequencies), np.fft.fftshift(amplitudes)


def find_mode_frequency(angular_frequencies, amplitudes,
                        approximate_mode_frequency, search_half_width):
    """Locate a solver mode and interpolate its physical positive frequency.

    The solver mode exp(-i omega_R t) appears at Fourier frequency -omega_R.
    We find the strongest power-spectrum bin near that approximate location,
    then fit a quadratic through log(power) at that bin and its two neighbours
    to estimate the peak between FFT bins.
    """
    angular_frequencies = np.asarray(angular_frequencies, dtype=float)
    amplitudes = np.asarray(amplitudes)
    if angular_frequencies.ndim != 1 or amplitudes.ndim != 1:
        raise ValueError(
            "angular_frequencies and amplitudes must both be one-dimensional")
    if angular_frequencies.size != amplitudes.size:
        raise ValueError(
            "angular_frequencies and amplitudes must have the same length")
    if not np.isfinite(approximate_mode_frequency):
        raise ValueError("approximate_mode_frequency must be finite")
    if not np.isfinite(search_half_width) or search_half_width <= 0.0:
        raise ValueError("search_half_width must be finite and positive")

    fourier_peak_guess = -approximate_mode_frequency
    search_region = np.abs(
        angular_frequencies - fourier_peak_guess) <= search_half_width
    candidate_indices = np.flatnonzero(search_region)
    if candidate_indices.size < 3:
        raise ValueError("frequency search region contains fewer than three bins")

    power = np.abs(amplitudes) ** 2
    peak_index = candidate_indices[np.argmax(power[candidate_indices])]
    if peak_index == 0 or peak_index == power.size - 1:
        raise ValueError("peak lies at the edge of the Fourier spectrum")

    local_log_power = np.log(np.maximum(
        power[peak_index - 1:peak_index + 2], np.finfo(float).tiny))
    curvature = (local_log_power[0] - 2.0 * local_log_power[1]
                 + local_log_power[2])
    if curvature >= 0.0:
        raise ValueError("selected spectral feature is not a local maximum")

    fractional_bin_offset = 0.5 * (
        local_log_power[0] - local_log_power[2]) / curvature
    frequency_bin_spacing = (
        angular_frequencies[peak_index + 1] - angular_frequencies[peak_index])
    interpolated_fourier_frequency = (
        angular_frequencies[peak_index]
        + fractional_bin_offset * frequency_bin_spacing)
    return -interpolated_fourier_frequency


def dolan_filter_profile(frequency_offset, filter_width):
    """Evaluate F(x) = exp[-x^4 / w^4] from Dolan Eq. (49).

    ``frequency_offset`` is omega - omega_hat for the mode being isolated, and
    ``filter_width`` is w.  Both use angular-frequency units.
    """
    if not np.isfinite(filter_width) or filter_width <= 0.0:
        raise ValueError("filter_width must be finite and positive")

    scaled_frequency_offset = np.asarray(frequency_offset) / filter_width
    return np.exp(-(scaled_frequency_offset ** 4))


def reconstruct_filtered_mode(times, signal, mode_angular_frequency,
                              filter_bin_count=5.0):
    """Reconstruct one solver mode using Dolan's frequency filter.

    ``mode_angular_frequency`` is the positive physical omega_R in the solver
    convention psi ~ exp(-i omega_R t).  Equation (46) therefore places its
    Fourier peak at -omega_R.  Following the text below Eq. (49), the default
    filter width is five Fourier bins: w = 5 Delta omega.
    """
    if not np.isfinite(mode_angular_frequency):
        raise ValueError("mode_angular_frequency must be finite")
    if not np.isfinite(filter_bin_count) or filter_bin_count <= 0.0:
        raise ValueError("filter_bin_count must be finite and positive")

    angular_frequencies, amplitudes = fourier_amplitudes(times, signal)
    sample_spacing = float(np.asarray(times)[1] - np.asarray(times)[0])
    frequency_bin_spacing = 2.0 * np.pi / (len(times) * sample_spacing)
    filter_width = filter_bin_count * frequency_bin_spacing

    frequency_offsets = angular_frequencies + mode_angular_frequency
    filter_values = dolan_filter_profile(frequency_offsets, filter_width)
    filtered_amplitudes = amplitudes * filter_values

    # fourier_amplitudes includes Dolan's 1/N normalization; NumPy's inverse
    # transform also divides by N, so restore that factor before inversion.
    unshifted_amplitudes = np.fft.ifftshift(filtered_amplitudes)
    return np.fft.ifft(unshifted_amplitudes * len(times))


def fit_growth_rate(times, filtered_mode, central_fraction=0.5):
    """Fit Gamma from the central part of a reconstructed mode.

    Dolan fits the mode square amplitude, for which
    log(|psi|^2) = constant + 2 Gamma t.  The ends of the reconstructed record
    depend on the frequency-filter shape, so the default fit keeps only the
    central half, centred on t_max/2 as recommended below Eq. (49).
    """
    times = np.asarray(times, dtype=float)
    filtered_mode = np.asarray(filtered_mode)
    if times.ndim != 1 or filtered_mode.ndim != 1:
        raise ValueError("times and filtered_mode must both be one-dimensional")
    if times.size != filtered_mode.size:
        raise ValueError("times and filtered_mode must have the same length")
    if not 0.0 < central_fraction <= 1.0:
        raise ValueError("central_fraction must lie in (0, 1]")

    record_midpoint = 0.5 * (times[0] + times[-1])
    half_fit_duration = 0.5 * central_fraction * (times[-1] - times[0])
    fit_region = np.abs(times - record_midpoint) <= half_fit_duration
    mode_square_amplitude = np.abs(filtered_mode[fit_region]) ** 2
    if np.any(mode_square_amplitude == 0.0):
        raise ValueError("filtered mode has zero amplitude inside the fit region")

    centered_fit_times = times[fit_region] - record_midpoint
    square_amplitude_slope, _ = np.polyfit(
        centered_fit_times, np.log(mode_square_amplitude), 1)
    return 0.5 * square_amplitude_slope
