"""Freeze the populated benchmark, reference branches, and requested domain."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tomllib
from freeze_source import freeze_source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    seeds = campaign / "work/direct_fig7_20261026_v5/inputs"
    sources = {
        "equal-energy.jls": campaign / "results/equal_energy_20261028_v1.jls",
        "primary.jls": seeds / "alpha030_cosine_branch.jls",
        "secondary.jls": seeds / "alpha030_cosine_secondary.jls",
    }
    for path in sources.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    with (campaign / "results/equal_energy_20261028_v1.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    freeze_source(args.name)
    snapshot = campaign / "work" / args.name
    inputs = snapshot / "inputs"
    inputs.mkdir()
    for name, path in sources.items():
        shutil.copy2(path, inputs / name)
    (inputs / "request.toml").write_text(
        f'reference_energy = {metadata["reference_energy"]!r}\n'
        'inner = -80.0\nouter = 180.0\nalpha = 0.3\nspin = 0.99\n'
        'joint_potentials = ["cosine", "chiral"]\nreference_potential = "cosine"\n'
        'physical_decay_constant_assigned = false\n')
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for path in inputs.iterdir():
        provenance["sha256"][str(path.relative_to(snapshot))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + name + "\n" for name, digest in sorted(provenance["sha256"].items())))
    print("Frozen files:", len(provenance["sha256"]))


if __name__ == "__main__":
    main()
