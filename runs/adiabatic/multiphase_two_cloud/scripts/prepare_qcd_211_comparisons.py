"""Freeze the remaining 23 first-peak, 211-dominated QCD comparisons."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
from freeze_source import freeze_source


def prepare(name):
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    population_path = campaign / "results/qcd_populations_20261029_v1/population_targets.csv"
    with population_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    cases, sources = [], {"population_targets.csv": population_path}
    for row in rows:
        alpha, spin, mass = (float(row[key]) for key in ("alpha", "spin", "mass_solar"))
        if int(row["dominant_m"]) != 1 or (alpha, spin, mass) == (.3, .99, 5.):
            continue
        label = f"qcd211_alpha{round(100*alpha):03d}_spin{round(100*spin):03d}_mass{round(mass):03d}"
        background = f"alpha{round(100*alpha):03d}_spin{round(100*spin):03d}"
        primary_name, secondary_name = f"{background}-primary.jls", f"{background}-secondary.jls"
        sources[primary_name] = (root / row["primary_reference"]).with_suffix(".jls")
        sources[secondary_name] = (root / row["secondary_reference"]).with_suffix(".jls")
        cases.append(dict(name=label, alpha=alpha, spin=spin, mass_solar=mass,
            dominant_m=1, dominant_mode="211", modes=[[1,1],[3,1],[5,1],[3,3],[5,3],[5,5]],
            polar_nodes=16, azimuthal_nodes=31, target_code_energy=float(row["dominant_code_energy"]),
            target_occupation=float(row["dominant_occupation"]),
            code_energy_to_mass_fraction=float(row["code_energy_to_mass_fraction"]),
            occupation_scale=float(row["occupation_scale"]),
            log10_minor_occupation=float(row["log10_minor_occupation"]),
            axion_ev=float(row["axion_ev"]), decay_constant_gev=float(row["decay_constant_gev"]),
            linear_seed=primary_name, secondary_seed=secondary_name,
            original_reference=row["primary_reference"]))
    assert len(cases) == 23
    assert all(source.is_file() for source in sources.values())
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
    (inputs / "cases.toml").write_text("\n".join(lines) + "\n")
    provenance = json.loads((snapshot / "provenance.json").read_text())
    provenance["sha256"].update({str(path.relative_to(snapshot)): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in inputs.iterdir()})
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        f"{digest}  {path}\n" for path, digest in sorted(provenance["sha256"].items())))
    print("Prepared", len(cases), "complete comparison jobs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    prepare(parser.parse_args().name)
