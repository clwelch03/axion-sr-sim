"""Screen whether the other cloud can grow before first-cloud GW depletion.

This evaluates algebraic rates and timescales, not a coupled time trajectory.
Use the mass-changing first-saturation endpoint, consistently expressing all
occupations and rates in units of the post-extraction BH mass. The chosen
cosine potential has lambda=mu^2/fa^2; c4=0.3 is retained as a comparison only.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil

from estimate_qcd_populations import (
    PLANCK_ENERGY_GEV, YEAR_SECONDS, saturation_endpoint, threshold_spin,
)

ROOT = Path(__file__).resolve().parents[4]
CAMPAIGN = ROOT / "runs/adiabatic/multiphase_two_cloud"


def occupation_rate(alpha, spin, azimuthal_number, *, binding=True):
    """Signed M*Gamma for n=m+1,l=m, with the hydrogenic horizon threshold.

    The factorial/product expression is the same leading small-alpha rate as
    common.hydrogenic. Optionally insert omega_R=mu*(1-alpha^2/(2*n^2)) in the
    horizon factors so the rate vanishes at the endpoint used in this screen.
    This replacement does not make the rate accurate at finite alpha.
    """
    degree = azimuthal_number
    principal = degree + 1
    frequency = alpha * (1-alpha**2/(2*principal**2)) if binding else alpha
    horizon = 1 + math.sqrt(1-spin**2)
    horizon_factor = azimuthal_number*spin - 2*frequency*horizon
    log_prefactor = ((4*degree+2)*math.log(2) + math.lgamma(2*degree+2)
                     - (2*degree+4)*math.log(principal)
                     + 2*(math.lgamma(degree+1)-math.lgamma(2*degree+1)
                          - math.lgamma(2*degree+2)))
    angular_product = math.prod(k*k*(1-spin*spin)+horizon_factor**2
                                for k in range(1, degree+1))
    return math.exp(log_prefactor)*alpha**(4*degree+5)*horizon_factor*angular_product


def rates(alpha, spin, decay_constant, quartic_coefficient=1.0):
    """M*gamma for epsilon=N/(G M^2); BGLS Tables IV--VII, PDF pp.49--50."""
    coupling = quartic_coefficient**2*(PLANCK_ENERGY_GEV/decay_constant)**4
    horizon = 1 + math.sqrt(1-spin**2)
    return dict(
        transfer_bh=4.3e-7*alpha**12*horizon*coupling,
        transfer_infinity=1.1e-8*alpha**9*coupling,
        annihilation_211=1e-2*alpha**15,
        annihilation_322=3e-8*alpha**19,
        transition_322_to_211=5e-6*alpha**11,
        relativistic_211=5e-9*alpha**22*coupling,
        relativistic_322=6e-14*alpha**28*coupling,
    )


def endpoint_coordinates(row, model):
    """Return alpha, spin, epsilon, N scale, and time unit at saturation."""
    alpha, spin = float(row["alpha"]), float(row["spin"])
    dominant = int(row["dominant_m"])
    if model == "mass_changing":
        endpoint = saturation_endpoint(alpha, spin, dominant)
        mass_ratio = endpoint["mass_ratio"]
        epsilon_initial_units = endpoint["epsilon"]
        endpoint_alpha, endpoint_spin = endpoint["alpha"], endpoint["spin"]
    elif model == "fixed_mass":
        mass_ratio = 1.0
        endpoint_alpha = alpha
        endpoint_spin = threshold_spin(alpha*(1-alpha**2/(2*(dominant+1)**2)), dominant)
        epsilon_initial_units = (spin-endpoint_spin)/dominant
    else:
        raise ValueError("Unknown endpoint model")
    return dict(alpha=endpoint_alpha, spin=endpoint_spin, mass_ratio=mass_ratio,
        epsilon=epsilon_initial_units/mass_ratio**2,
        occupation_scale=float(row["occupation_scale"])*mass_ratio**2,
        time_unit_seconds=float(row["time_unit_seconds"])*mass_ratio,
        epsilon_initial_units=epsilon_initial_units)


def log_minor_seed_at_saturation(row, epsilon_initial_units, *, seed_occupation=1.0,
                                  quartic_coefficient=1.0, transfer_multiplier=1.0):
    """Use the earlier growing-phase diagnostic only to set the logarithm.

    A damped other mode keeps a microscopic reference seed. This is not a
    prediction of a precise microscopic population. Vary that reference below.
    """
    dominant = int(row["dominant_m"])
    dominant_rate = float(row[f"occupation_rate_{'211' if dominant == 1 else '322'}"])
    minor_rate = float(row[f"occupation_rate_{'322' if dominant == 1 else '211'}"])
    if minor_rate <= 0:
        return math.log(seed_occupation)
    peak_occupation = epsilon_initial_units*float(row["occupation_scale"])
    growth_duration = math.log(peak_occupation/seed_occupation)/dominant_rate
    initial_rates = rates(float(row["alpha"]), float(row["spin"]),
                          float(row["decay_constant_gev"]), quartic_coefficient)
    pump = initial_rates["transfer_bh" if dominant == 1 else "transfer_infinity"]
    seed_epsilon = seed_occupation/float(row["occupation_scale"])
    exponent = transfer_multiplier*pump*(epsilon_initial_units**2-seed_epsilon**2)/(2*dominant_rate)
    return math.log(seed_occupation)+minor_rate*growth_duration+exponent


def growth_exponent_with_gw_depletion(time, half_life, horizon_rate, pump_rate,
                                     transition_gain=0.0):
    """Analytic gain for epsilon_dom(t)=epsilon_0/(1+t/tau_GW).

    Pump falls as epsilon_dom^2; a positive stimulated GW transition gain
    falls as epsilon_dom. This is a linearized, fixed-background diagnostic,
    not a coupled population evolution or a prediction after catch-up.
    """
    ratio = time/half_life
    return (horizon_rate*time + pump_rate*half_life*ratio/(1+ratio)
            + transition_gain*half_life*math.log1p(ratio))


def screen(row, *, endpoint_model="mass_changing", quartic_coefficient=1.0,
           relative_growth_multiplier=1.0, seed_occupation=1.0,
           catchup_fraction=1.0, binding=True):
    coordinates = endpoint_coordinates(row, endpoint_model)
    alpha, spin, epsilon = (coordinates[key] for key in ("alpha", "spin", "epsilon"))
    dominant = int(row["dominant_m"])
    minor = 3-dominant
    coefficients = rates(alpha, spin, float(row["decay_constant_gev"]), quartic_coefficient)
    horizon_rate = occupation_rate(alpha, spin, minor, binding=binding)
    pump_rate = coefficients["transfer_bh" if dominant == 1 else "transfer_infinity"]*epsilon**2
    # Keep gains that could overturn a damping exclusion. Omitting the reverse
    # transition loss for a 322 minor mode makes the screening test optimistic.
    transition_gain = coefficients["transition_322_to_211"]*epsilon if dominant == 2 else 0.0
    horizon_rate *= relative_growth_multiplier
    pump_rate *= relative_growth_multiplier
    transition_gain *= relative_growth_multiplier
    net_rate = horizon_rate+pump_rate+transition_gain
    annihilation = coefficients[f"annihilation_{'211' if dominant == 1 else '322'}"]
    half_life = 1/(2*annihilation*epsilon)
    peak_occupation = epsilon*coordinates["occupation_scale"]
    log_seed = log_minor_seed_at_saturation(row, coordinates["epsilon_initial_units"],
        seed_occupation=seed_occupation, quartic_coefficient=quartic_coefficient)
    required_efolds = math.log(catchup_fraction*peak_occupation)-log_seed
    if required_efolds <= 0:
        raise ValueError("Seed is already macroscopic; this screen does not apply")
    margin = net_rate*half_life/required_efolds
    growth_time = required_efolds/net_rate if net_rate > 0 else None
    years = coordinates["time_unit_seconds"]/YEAR_SECONDS
    if net_rate <= 0:
        status = "minor_net_damped"
        half_life_gain = 0.0
    else:
        status = "passes_optimistic_timescale_test" if margin >= 1 else "growth_slower_than_depletion"
        half_life_gain = max(0.0, growth_exponent_with_gw_depletion(
            half_life, half_life, horizon_rate, pump_rate, transition_gain))
    # The comparison threshold at one half-life is fraction*N_dom(0)/2.
    half_life_required = required_efolds-math.log(2)
    scalar_loss = 3*coefficients[f"relativistic_{'211' if dominant == 1 else '322'}"]*epsilon**2
    return dict(
        status=status, passes=net_rate > 0 and margin >= 1,
        endpoint_model=endpoint_model, endpoint_alpha=alpha, endpoint_spin=spin,
        endpoint_mass_ratio=coordinates["mass_ratio"],
        epsilon_in_endpoint_mass_units=epsilon, dominant_occupation=peak_occupation,
        log_minor_seed_at_saturation=log_seed, required_efolds=required_efolds,
        minor_horizon_rate=horizon_rate, minor_pump_rate=pump_rate,
        minor_positive_gw_transition_rate=transition_gain, minor_net_rate=net_rate,
        gw_annihilation_coefficient=annihilation,
        gw_depletion_half_life_years=half_life*years,
        optimistic_catchup_years=growth_time*years if growth_time is not None else None,
        depletion_over_catchup_ratio=max(0.0, margin),
        relative_growth_enhancement_to_pass=1/margin if margin > 0 else None,
        integrated_gain_fraction_at_half_life=half_life_gain/half_life_required,
        passes_analytic_half_life_test=half_life_gain >= half_life_required,
        scalar_emission_loss_over_gw_loss=scalar_loss*half_life,
        pump_and_transition_over_horizon_damping=(pump_rate+transition_gain)/abs(horizon_rate)
            if horizon_rate < 0 else None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    population_path = CAMPAIGN / "results/qcd_populations_20261029_v1/population_targets.csv"
    with population_path.open(newline="") as stream:
        inputs = list(csv.DictReader(stream))
    assert len(inputs) == 36
    rows, sensitivity = [], []
    for source in inputs:
        identity = {key:float(source[key]) for key in ("alpha", "spin", "mass_solar")}
        identity.update(dominant_mode=source["dominant_mode"], minor_mode=source["minor_mode"])
        nominal = screen(source)
        endpoint = saturation_endpoint(identity["alpha"], identity["spin"], int(source["dominant_m"]))
        assert math.isclose(endpoint["epsilon"], float(source["endpoint_normalized_occupation"]), rel_tol=1e-12)
        variations = []
        for model, seed, fraction, binding in itertools.product(
                ("mass_changing", "fixed_mass"), (1.0, 1e6), (1.0, 0.1), (True, False)):
            result = screen(source, endpoint_model=model, seed_occupation=seed,
                            catchup_fraction=fraction, binding=binding)
            variations.append(result)
            sensitivity.append(dict(**identity, seed_occupation=seed,
                catchup_fraction=fraction, binding_horizon_factor=binding, **result))
        qcd = screen(source, quartic_coefficient=0.3)
        enhanced = screen(source, relative_growth_multiplier=10.0)
        rows.append(dict(**identity, **nominal,
            alternative_assumptions_any_pass=any(item["passes"] for item in variations),
            alternative_assumptions_max_ratio=max(item["depletion_over_catchup_ratio"] for item in variations),
            qcd_c4_03_ratio=qcd["depletion_over_catchup_ratio"],
            tenfold_growth_over_gw_pass=enhanced["passes"],
            tenfold_growth_analytic_half_life_pass=enhanced["passes_analytic_half_life_test"]))
    args.output_directory.mkdir(parents=True, exist_ok=False)
    for name, records in (("catchup_screen.csv", rows), ("assumption_sensitivity.csv", sensitivity)):
        with (args.output_directory/name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    dependencies = [Path(__file__), Path(__file__).with_name("estimate_qcd_populations.py"),
                    ROOT / "src/common/hydrogenic.py", population_path]
    metadata = dict(
        nominal_passes=sum(row["passes"] for row in rows), cases=len(rows),
        damped_cases=sum(row["status"] == "minor_net_damped" for row in rows),
        tenfold_growth_sensitivity_passes=sum(row["tenfold_growth_over_gw_pass"] for row in rows),
        potential="cosine", quartic_coefficient=1.0, effective_quartic_scale="f=fa",
        qcd_relation="mu[eV]*fa[GeV]=5.7e6", age_cut_imposed=False,
        criterion="ln(N_dom_peak/N_minor_at_saturation)/(Gamma_minor+gamma_pump*epsilon_dom^2+positive_GW_transition) <= 1/(2*gamma_GW_dom*epsilon_dom)",
        rate_units="all tabulated rates dimensionless in post-extraction BH time units G*M_final/c^3",
        endpoint="algebraic mass-changing first saturation, hydrogenic real frequency; epsilon=N/(G*M_final^2)",
        definitions=dict(catchup="minor reaches initial dominant occupation in optimistic frozen-rate estimate",
            depletion="GW-only half-life of the dominant cloud; additional losses shorten it",
            sensitivity="tenfold faster net small-cloud growth relative to GW depletion, not a confidence interval"),
        sources=[dict(url="https://arxiv.org/pdf/2011.11646v2",
            references="PDF Eqs. (39), (42)-(43); Tables IV-VII; Appendix F",
            role="GW depletion criterion, leading rates, and saturation endpoint")],
        limitations=[
            "Necessary optimistic timescale screen, not a demonstration of sustained coexistence or a coupled time trajectory.",
            "Pump held at its initial maximum in nominal timescale test; analytic half-life column lets it decline under GW-only depletion.",
            "Negative horizon rates are retained. Losses of a 322 minor mode through GW transitions/mixed processes are omitted, making catch-up easier.",
            "Binding correction only adjusts the horizon factor; no numerical post-extraction eigenfrequencies were computed.",
            "Small-alpha interaction/GW rates and a two-level truncation are preliminary at alpha=0.2-0.4. Close failures require more accurate rates.",
            "No BH lifetime is imposed. A later dominant cloud after the original has depleted is not excluded.",
            "Original profile population targets and submitted jobs are unchanged by this diagnostic.",
        ],
        sha256={str(path.resolve().relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in dependencies})
    (args.output_directory/"metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    for source in dependencies[:2]:
        shutil.copy2(source, args.output_directory/source.name)
    print(json.dumps({key:metadata[key] for key in ("cases", "nominal_passes", "damped_cases",
        "tenfold_growth_sensitivity_passes")}, indent=2))
    for row in rows:
        if row["mass_solar"] == 5:
            print(row["alpha"], row["spin"], row["dominant_mode"], row["status"],
                  "depletion/catchup=",row["depletion_over_catchup_ratio"])


if __name__ == "__main__":
    main()
