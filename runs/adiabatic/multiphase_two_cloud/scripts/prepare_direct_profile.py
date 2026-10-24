"""Freeze source and the three saved Omiya inputs for one direct profile solve."""

import hashlib
import json
from pathlib import Path
import shutil
import sys

from freeze_source import freeze_source


def prepare(name):
    root = Path(__file__).resolve().parents[4]
    reference = root / "runs/adiabatic/multimode_secondary/results"
    names = (
        "alpha030_cosine_branch.jls",
        "alpha030_cosine_secondary.jls",
        "alpha030_cosine_secondary_energy.csv",
    )
    for filename in names:
        if not (reference / filename).is_file():
            raise FileNotFoundError(reference / filename)
    freeze_source(name)
    snapshot = root / "runs/adiabatic/multiphase_two_cloud/work" / name
    inputs = snapshot / "inputs"
    inputs.mkdir()
    manifest_path = snapshot / "provenance.json"
    manifest = json.loads(manifest_path.read_text())
    for filename in names:
        destination = inputs / filename
        shutil.copy2(reference / filename, destination)
        manifest["sha256"][str(destination.relative_to(snapshot))] = hashlib.sha256(
            destination.read_bytes()
        ).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("Source and reference inputs frozen. No job submitted.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_direct_profile.py NEW_RUN_NAME")
    prepare(sys.argv[1])
