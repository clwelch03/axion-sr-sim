"""Freeze current src and the verified inputs for one target collocation run."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
from freeze_source import freeze_source

def prepare(name):
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    inputs = {
        "joint-input.jls": campaign / "results/direct_fig7_20261026_v5-input.jls",
        "initial.jls": root / "runs/adiabatic/radial_bvp/results/target_driven_mirk4_refined_v1.jls",
        "reference-curves.jls": campaign / "results/reference_curves_20261028_v1.jls",
        "reference-curves.csv": campaign / "results/reference_curves_20261028_v1.csv",
        "reference-curves.toml": campaign / "results/reference_curves_20261028_v1.toml",
    }
    for source in inputs.values():
        if not source.is_file():
            raise FileNotFoundError(source)
    freeze_source(name)
    snapshot = campaign / "work" / name
    (snapshot / "inputs").mkdir()
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for filename, source in inputs.items():
        target = snapshot / "inputs" / filename
        shutil.copy2(source, target)
        provenance["sha256"][str(target.relative_to(snapshot))] = hashlib.sha256(target.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + name + "\n" for name, digest in sorted(provenance["sha256"].items())))
    print("Prepared target snapshot:", snapshot)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_collocation_target.py UNIQUE_SNAPSHOT_NAME")
    prepare(sys.argv[1])
