"""Estimate first-saturation QCD populations without time integration.

Use the saved linear frequencies, the leading spin-extraction occupation,
and the perturbative transfer coefficients of Baryakhtar et al. (2021).
The second cloud is seeded at occupation one for a growing-phase diagnostic;
microscopic occupations are not claimed as accurate classical predictions.
This script prepares a catalogue, not nonlinear solver submissions.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tomllib

import numpy as np
from scipy.constants import G, c, hbar, electron_volt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from common.hydrogenic import growth_rate

SOLAR_MASS_KG = 1.988409870698051e30
PLANCK_MASS_KG = math.sqrt(hbar * c / G)  # unreduced
PLANCK_ENERGY_GEV = PLANCK_MASS_KG * c**2 / electron_volt / 1e9
YEAR_SECONDS = 365.25 * 86400
QCD_QUARTIC_COEFFICIENT = 0.3  # lambda = 0.3 * mu^2 / fa^2
BACKGROUND_GRID = tuple(itertools.product((0.2, 0.3, 0.4), (0.8, 0.9, 0.95, 0.99)))
JULIA = "/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia"


def threshold_spin(frequency_over_mass, azimuthal_number):
    """Invert omega = m Omega_H; the argument is dimensionless M*omega."""
    ratio = frequency_over_mass / azimuthal_number
    return 4 * ratio / (1 + 4 * ratio**2)


def leading_saturation(alpha, spin, azimuthal_number):
    """Return epsilon=N/(G M^2), with hbar=c=1, at leading spin extraction."""
    return max(0.0, spin - threshold_spin(alpha, azimuthal_number)) / azimuthal_number


def saturation_endpoint(alpha, spin, azimuthal_number):
    """Mass-changing hydrogenic endpoint, generalized from BGLS Appendix F.

    Eliminate time from dM=-omega*dN and dJ=-m*dN. For n=m+1,
    omega/mu=1-alpha_final^2/(2*n^2). Angular-momentum conservation then
    gives spin_final=(spin-m*epsilon)/(M_final/M_initial)^2. A scalar root
    locates the SR threshold. This is a sensitivity estimate, not a new target
    background or a calculation including cloud self-gravity and radiation.
    """
    principal_number = azimuthal_number + 1
    binding_frequency = alpha * (1 - alpha**2 / (2 * principal_number**2))
    if spin <= threshold_spin(binding_frequency, azimuthal_number):
        return dict(epsilon=0.0, mass_ratio=1.0, alpha=alpha, spin=spin,
                    cloud_mass_fraction=0.0)
    scale = math.sqrt(2) * principal_number

    def final_alpha(epsilon):
        return scale * math.tanh(math.atanh(alpha / scale) - alpha**2 * epsilon / scale)

    def final_spin(epsilon):
        mass_ratio = final_alpha(epsilon) / alpha
        return (spin - azimuthal_number * epsilon) / mass_ratio**2

    def residual(epsilon):
        background = final_alpha(epsilon)
        frequency = background * (1 - background**2 / (2 * principal_number**2))
        return final_spin(epsilon) - threshold_spin(frequency, azimuthal_number)

    epsilon = brentq(residual, 0.0, spin / azimuthal_number, xtol=1e-14)
    background = final_alpha(epsilon)
    return dict(epsilon=epsilon, mass_ratio=background / alpha, alpha=background,
                spin=final_spin(epsilon), cloud_mass_fraction=epsilon * alpha *
                (1 - background**2 / (2 * principal_number**2)))


def physical_units(alpha, mass_solar):
    mass = mass_solar * SOLAR_MASS_KG
    time_unit = G * mass / c**3
    axion_ev = alpha * hbar / time_unit / electron_volt
    decay_constant = 5.7e6 / axion_ev  # mu[eV] * fa[GeV] = 5.7e6
    return dict(time_unit_seconds=time_unit, axion_ev=axion_ev,
                decay_constant_gev=decay_constant,
                occupation_scale=(mass / PLANCK_MASS_KG)**2,
                code_energy_to_mass_fraction=(decay_constant / PLANCK_ENERGY_GEV)**2)


def transfer_rates(alpha, spin, decay_constant, quartic_coefficient=QCD_QUARTIC_COEFFICIENT):
    """Return M*gamma for epsilon occupations; Table I and Table VI of BGLS.

    Their f is the quartic scale, f=fa/sqrt(c4), not the QCD decay constant.
    Converting gamma/mu to M*gamma supplies one extra power of alpha.
    """
    coupling = quartic_coefficient**2 * (PLANCK_ENERGY_GEV / decay_constant)**4
    into_bh = 4e-7 * alpha**12 * (1 + math.sqrt(1 - spin**2)) * coupling
    to_infinity = 1.1e-8 * alpha**9 * coupling
    return into_bh, to_infinity


def minor_growth_diagnostic(epsilon, occupation_scale, dominant_rate, minor_rate,
                            pump_rate, seed_occupation=1.0):
    """Analytic linearized growth with epsilon_dom exponentially increasing.

    log(N_minor/N_seed) = gamma_minor*t + gamma_pump/(2*gamma_dom)
    * (epsilon_dom^2-epsilon_seed^2). This extends BGLS Eq. (E3) by retaining
    independent linear growth. For a damped minor mode, report the peak gain
    versus damping; do not extrapolate sub-particle occupations as predictions.
    """
    if dominant_rate <= 0 or epsilon * occupation_scale <= seed_occupation:
        raise ValueError("Require a growing dominant cloud above its seed")
    e_folds = math.log(epsilon * occupation_scale / seed_occupation)
    dimensionless_time = e_folds / dominant_rate
    transfer_exponent = pump_rate * (epsilon**2 - (seed_occupation / occupation_scale)**2) / (
        2 * dominant_rate)
    peak_minor_rate = minor_rate + pump_rate * epsilon**2
    if minor_rate <= 0 and peak_minor_rate > 0:
        # A damped mode becomes unstable: the vacuum-seeded early-time formula
        # cannot determine its occupation once the sign changes.
        return dict(log10_minor_occupation=None, log10_minor_over_dominant=None,
                    dimensionless_time=dimensionless_time, transfer_exponent=transfer_exponent,
                    peak_minor_rate=peak_minor_rate, status="minor_instability_unresolved")
    if minor_rate <= 0:
        log_minor = math.log(seed_occupation)  # microscopic reference, not an imposed floor in an ODE
        status = "minor_damped_microscopic_reference"
    else:
        log_minor = math.log(seed_occupation) + minor_rate * dimensionless_time + transfer_exponent
        status = "growing_phase_diagnostic"
    return dict(log10_minor_occupation=log_minor / math.log(10),
                log10_minor_over_dominant=(log_minor-math.log(epsilon*occupation_scale))/math.log(10),
                dimensionless_time=dimensionless_time, transfer_exponent=transfer_exponent,
                peak_minor_rate=peak_minor_rate, status=status)


def load_linear_references(results_directory, julia):
    """Read 23 TOML references and the separately recovered spin-0.95 JLS."""
    references, inputs = {}, set()
    for alpha, spin in BACKGROUND_GRID:
        name = f"alpha{round(100*alpha):03d}_spin{round(100*spin):03d}"
        for cloud, azimuthal_number in (("primary", 1), ("secondary", 2)):
            filename = name + "-" + cloud + ".toml"
            candidates = [results_directory / directory / filename for directory in
                          ("linear_grid_retry_20261029_v1", "linear_grid_20261029_v1")]
            path = next((path for path in candidates if path.is_file()), None)
            if path is None:
                if (alpha, spin, azimuthal_number) != (0.3, 0.95, 2):
                    raise FileNotFoundError(filename)
                directory = results_directory / "omiya_spin095_linear_20261029_v2"
                path = directory / "linear_reference.jls"
                metadata_path = path.with_suffix(".toml")
                metadata = tomllib.loads(metadata_path.read_text())
                if not metadata["passed"]:
                    raise ValueError("Recovered reference failed")
                code = ('using Serialization, TOML; x=deserialize(ARGS[1]); '
                        's=x.references["linear_secondary"].solution; '
                        'TOML.print(stdout, Dict("frequency_real"=>real(s.frequency), '
                        '"frequency_imaginary"=>imag(s.frequency)))')
                values = tomllib.loads(subprocess.check_output(
                    [julia, "--startup-file=no", "--compiled-modules=existing", "-e", code, str(path)],
                    text=True))
                factor = metadata["amplitude_factors"]["linear_secondary"]
                values["energy"] = metadata["reference_energy"] / factor**2
                inputs.add(metadata_path)
            else:
                values = tomllib.loads(path.read_text())
                if not values["passed"]:
                    raise ValueError("Reference failed: " + str(path))
            inputs.add(path)
            references[alpha, spin, azimuthal_number] = dict(
                frequency=values["frequency_real"], rate=2*values["frequency_imaginary"],
                unit_outer_energy=values["energy"], source=str(path.relative_to(ROOT)))
    return references, inputs


def estimate_case(alpha, spin, mass_solar, references):
    modes = [references[alpha, spin, m] for m in (1, 2)]
    dominant = 1 + max(range(2), key=lambda index: modes[index]["rate"])
    minor = 3 - dominant
    unit = physical_units(alpha, mass_solar)
    epsilon = leading_saturation(alpha, spin, dominant)
    endpoint = saturation_endpoint(alpha, spin, dominant)
    if epsilon <= 0:
        raise ValueError("Leading threshold disagrees with selected growing mode")
    bh_pump, infinity_pump = transfer_rates(alpha, spin, unit["decay_constant_gev"])
    pump = bh_pump if dominant == 1 else infinity_pump
    analytic_rates = [2*float(growth_rate(alpha, spin, l=m, m=m)) for m in (1, 2)]
    numerical_rates = [mode["rate"] for mode in modes]
    if 1 + int(np.argmax(analytic_rates)) != dominant:
        raise ValueError("Analytic and saved linear modes disagree on dominant cloud")
    diagnostic = minor_growth_diagnostic(epsilon, unit["occupation_scale"],
        numerical_rates[dominant-1], numerical_rates[minor-1], pump)
    variations = []
    # This is a sensitivity sample, not a confidence interval: use both rate
    # estimates, both population estimates, 0.1x/10x transfer rates, and seeds
    # of one or a million particles (both negligible compared with ~1e77).
    for rates, saturation, factor, seed in itertools.product(
            (numerical_rates, analytic_rates), (epsilon, endpoint["epsilon"]),
            (0.1, 1.0, 10.0), (1.0, 1e6)):
        variations.append(minor_growth_diagnostic(saturation, unit["occupation_scale"],
            rates[dominant-1], rates[minor-1], factor*pump, seed))
    unresolved = any(item["log10_minor_over_dominant"] is None for item in variations)
    maximum_ratio = None if unresolved else max(item["log10_minor_over_dominant"] for item in variations)
    negligible = maximum_ratio is not None and maximum_ratio < -12
    mass_fraction = epsilon * modes[dominant-1]["frequency"]
    code_energy = mass_fraction / unit["code_energy_to_mass_fraction"]
    target_energies = [0.0, 0.0]
    target_energies[dominant-1] = code_energy
    if not negligible:
        target_energies[minor-1] = None
    row = dict(alpha=alpha, spin=spin, mass_solar=mass_solar,
        dominant_mode="211" if dominant == 1 else "322", dominant_m=dominant,
        minor_mode="322" if dominant == 1 else "211", **unit,
        dominant_normalized_occupation=epsilon,
        dominant_occupation=epsilon*unit["occupation_scale"],
        dominant_mass_fraction=mass_fraction, dominant_code_energy=code_energy,
        target_energy_211=target_energies[0], target_energy_322=target_energies[1],
        dominant_outer_amplitude_seed=math.sqrt(code_energy/modes[dominant-1]["unit_outer_energy"]),
        endpoint_normalized_occupation=endpoint["epsilon"],
        endpoint_over_leading_occupation=endpoint["epsilon"]/epsilon,
        endpoint_mass_fraction=endpoint["cloud_mass_fraction"],
        endpoint_alpha=endpoint["alpha"], endpoint_spin=endpoint["spin"],
        endpoint_mass_ratio=endpoint["mass_ratio"],
        occupation_rate_211=numerical_rates[0], occupation_rate_322=numerical_rates[1],
        hydrogenic_rate_211=analytic_rates[0], hydrogenic_rate_322=analytic_rates[1],
        transfer_rate_into_bh=bh_pump, transfer_rate_to_infinity=infinity_pump,
        **diagnostic,
        fixed_rate_growth_years=diagnostic["dimensionless_time"]*unit["time_unit_seconds"]/YEAR_SECONDS,
        sensitivity_max_log10_minor_over_dominant=maximum_ratio,
        minor_classical_population_negligible=negligible,
        minor_population_within_factor_ten_established=False,
        recommendation="single_dominant_cloud_plus_linear_probe" if negligible else "population_review_needed",
        primary_reference=modes[0]["source"], secondary_reference=modes[1]["source"])
    # At a 211 peak this is the competition with GW annihilation during growth.
    row["primary_gw_loss_over_growth"] = (
        2e-2*alpha**15*epsilon/numerical_rates[0] if dominant == 1 else None)
    row["minor_pump_over_damping_at_endpoint"] = (
        pump*endpoint["epsilon"]**2/abs(numerical_rates[minor-1])
        if numerical_rates[minor-1] < 0 else None)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path, help="new directory; existing outputs are preserved")
    parser.add_argument("--julia", default=JULIA, help="Julia executable for the one recovered JLS reference")
    args = parser.parse_args()
    results = ROOT / "runs/adiabatic/multiphase_two_cloud/results"
    references, inputs = load_linear_references(results, args.julia)
    rows = [estimate_case(alpha, spin, mass, references)
            for alpha, spin in BACKGROUND_GRID for mass in (5.0, 20.0, 100.0)]
    if len(rows) != 36 or any(not math.isfinite(row["dominant_code_energy"]) for row in rows):
        raise ValueError("Invalid catalogue")
    args.output_directory.mkdir(parents=True, exist_ok=False)
    with (args.output_directory / "population_targets.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    inputs.update((Path(__file__).resolve(), ROOT / "src/common/hydrogenic.py"))
    metadata = dict(
        stage="first dominant saturation scale, no pre-existing classical cloud",
        background="original grid retained; endpoint columns are sensitivity estimates only",
        rate_time_coordinate="t/(G*M/c^3); all rates are occupation rates, 2*omega_I",
        seed_occupation=1.0, qcd_quartic_coefficient=QCD_QUARTIC_COEFFICIENT,
        qcd_relation="mu[eV] * fa[GeV] = 5.7e6; quartic f = fa/sqrt(0.3)",
        minor_zero_convention="zero target is the negligible-classical-cloud limit, not literal N=0",
        sources=[
            dict(url="https://arxiv.org/pdf/1411.2263v3", equations="(8)", role="leading spin budget"),
            dict(url="https://arxiv.org/pdf/2011.11646v2", equations="(32), (34)-(35), (E1)-(E3), Appendix F, Tables I and VI", role="transfer rates and analytic estimates"),
            dict(url="https://arxiv.org/pdf/1511.02867", role="QCD mass-decay-constant relation"),
        ],
        limitations=[
            "Finite-population profile solves, time trajectories, and QCD-potential implementation were not performed.",
            "Linear references do not establish tiny-growth-rate or spatial convergence.",
            "Transfer coefficients are leading small-alpha approximations; higher modes, self-gravity and late coexistence are not modeled.",
            "The early-time exponent is only a diagnostic at saturation; mass/spin drift and GW depletion change the late growth.",
            "Seed dependence precludes a factor-ten prediction for microscopic minor occupations; sensitivity ranges are not error bars.",
            "Single-cloud target energy uses the saved linear real frequency on the original background.",
        ],
        sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(inputs)},
        python=sys.version, numpy=np.__version__)
    (args.output_directory / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    (args.output_directory / "estimate_qcd_populations.py").write_bytes(Path(__file__).read_bytes())
    print(f"Saved {len(rows)} cases to {args.output_directory}")
    for alpha, spin in BACKGROUND_GRID:
        row = next(row for row in rows if row["alpha"]==alpha and row["spin"]==spin)
        print(f"alpha={alpha} spin={spin}: {row['dominant_mode']}, "
              f"mass fraction={row['dominant_mass_fraction']:.5g}, "
              f"endpoint/leading={row['endpoint_over_leading_occupation']:.3g}")
    print("Negligible minor cloud:", sum(row["minor_classical_population_negligible"] for row in rows), "/ 36")


if __name__ == "__main__":
    main()
