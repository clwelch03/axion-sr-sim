"""Freeze a joint-profile sweep with explicit, independent energy-pair inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tomllib

from freeze_source import freeze_source


def prepare(name, energy_pairs):
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    seed = campaign / "results/equal_energy_20261028_v1.jls"
    metadata = campaign / "results/equal_energy_20261028_v1.toml"
    baseline = tomllib.loads(metadata.read_text())
    if not baseline["passed"] or not seed.is_file():
        raise ValueError("A converged equal-energy seed is required")
    if not energy_pairs or any(len(pair) != 2 or any(not 0 < value < float("inf") for value in pair)
                               for pair in energy_pairs):
        raise ValueError("Supply positive finite (E1/E0, E2/E0) pairs")
    if len(set(map(tuple, energy_pairs))) != len(energy_pairs):
        raise ValueError("Duplicate energy pairs")
    freeze_source(name, seed)
    snapshot = campaign / "work" / name
    shutil.copy2(metadata, snapshot / "inputs/initial.toml")
    reference_energy = baseline["reference_energy"]
    lines = [f"reference_energy = {reference_energy!r}",
             'energy_units = "F_a^2 M"',
             'normalization = "independent_cloud_energy_targets"',
             "# Domains are provisional; endpoint/angular convergence remains separate.", ""]
    case_count = 0
    # Start with the existing background, then move outward in parameter space.
    for alpha in (0.3, 0.2, 0.4):
        for spin in (0.99, 0.95, 0.9, 0.8):
            outer = 675.0 if alpha == 0.2 else 300.0
            for pair_index, pair in enumerate(energy_pairs, 1):
                case_count += 1
                case_name = f"alpha{round(100*alpha):03d}_spin{round(100*spin):03d}_energy{pair_index}"
                energies = [factor * reference_energy for factor in pair]
                lines.extend(["[[cases]]", f'name = "{case_name}"', f"alpha = {alpha!r}",
                              f"spin = {spin!r}", f"energy_factors = {list(pair)!r}",
                              f"energies = {energies!r}", f"outer = {outer!r}", ""])
    config_path = snapshot / "inputs/sweep.toml"
    config_path.write_text("\n".join(lines))
    parsed = tomllib.loads(config_path.read_text())
    assert len(parsed["cases"]) == case_count
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for path in sorted((snapshot / "inputs").iterdir()):
        provenance["sha256"][str(path.relative_to(snapshot))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + path + "\n" for path, digest in sorted(provenance["sha256"].items())))
    print(f"Prepared {case_count} configurations: {config_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("--energy-pair", type=float, nargs=2, action="append", required=True,
                        metavar=("E1_OVER_E0", "E2_OVER_E0"))
    args = parser.parse_args()
    prepare(args.name, args.energy_pair)
