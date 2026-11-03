"""Freeze the two QCD pilot cases for direct, fixed-amplitude nonlinear solves."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tomllib

from freeze_source import freeze_source


def prepare(name):
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    population_directory = campaign / "results/qcd_populations_20261029_v1"
    population_path = population_directory / "population_targets.csv"
    with population_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    cases, sources = [], {}
    for alpha, spin, dominant_m in ((0.3, 0.99, 1), (0.4, 0.95, 2)):
        matches = [row for row in rows if float(row["alpha"]) == alpha
                   and float(row["spin"]) == spin and float(row["mass_solar"]) == 5]
        if len(matches) != 1:
            raise ValueError("Expected exactly one population row per pilot")
        row = matches[0]
        assert int(row["dominant_m"]) == dominant_m
        label = f"qcd{row['dominant_mode']}_alpha{round(100*alpha):03d}_spin{round(100*spin):03d}_mass005"
        reference = root / row["primary_reference" if dominant_m == 1 else "secondary_reference"]
        with reference.open("rb") as stream:
            report = tomllib.load(stream)
        assert report["passed"] and report["m"] == dominant_m
        assert report["normalization"] == "first outer coefficient = 1"
        sources[label + "-linear.jls"] = reference.with_suffix(".jls")
        sources[label + "-linear.toml"] = reference
        modes = ([[1, 1], [3, 1], [5, 1], [3, 3], [5, 3], [5, 5]] if dominant_m == 1
                 else [[2, 2], [4, 2], [6, 2], [6, 6], [8, 6], [10, 10]])
        cases.append(dict(name=label, alpha=alpha, spin=spin, mass_solar=5.0,
            dominant_m=dominant_m, dominant_mode=row["dominant_mode"], modes=modes,
            polar_nodes=16 if dominant_m == 1 else 24,
            azimuthal_nodes=31 if dominant_m == 1 else 61,
            target_code_energy=float(row["dominant_code_energy"]),
            target_occupation=float(row["dominant_occupation"]),
            code_energy_to_mass_fraction=float(row["code_energy_to_mass_fraction"]),
            axion_ev=float(row["axion_ev"]), decay_constant_gev=float(row["decay_constant_gev"]),
            linear_seed=label + "-linear.jls", original_reference=str(reference.relative_to(root))))
    sources["population_targets.csv"] = population_path
    sources["population_metadata.json"] = population_directory / "metadata.json"
    for source in sources.values():
        if not source.is_file():
            raise FileNotFoundError(source)
    freeze_source(name)
    snapshot = campaign / "work" / name
    inputs = snapshot / "inputs"
    inputs.mkdir()
    for destination, source in sources.items():
        shutil.copy2(source, inputs / destination)
    lines = ["matching_tolerance = 1e-8", "max_iterations = 15", "relative_energy_window = 0.1"]
    for case in cases:
        lines.append("\n[[cases]]")
        lines.extend(f"{key} = {json.dumps(value)}" for key, value in case.items())
    (inputs / "pilots.toml").write_text("\n".join(lines) + "\n")
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for path in inputs.iterdir():
        provenance["sha256"][str(path.relative_to(snapshot))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + path + "\n" for path, digest in sorted(provenance["sha256"].items())))
    print("Prepared two direct QCD pilot solves:", snapshot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    prepare(parser.parse_args().name)
