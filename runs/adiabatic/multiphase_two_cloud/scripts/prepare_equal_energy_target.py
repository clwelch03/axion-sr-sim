"""Freeze source, the converged joint solution, and equal-energy references."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
from freeze_source import freeze_source


def prepare(name):
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    sources = {
        "joint-input.jls": campaign / "results/direct_fig7_20261026_v5-input.jls",
        "initial.jls": campaign / "results/src_target_20261028_v2.jls",
        "reference-curves.jls": campaign / "results/reference_curves_20261028_v2.jls",
        "reference-curves.toml": campaign / "results/reference_curves_20261028_v2.toml",
    }
    for path in sources.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    freeze_source(name)
    snapshot = campaign / "work" / name
    (snapshot / "inputs").mkdir()
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for name, path in sources.items():
        target = snapshot / "inputs" / name
        shutil.copy2(path, target)
        provenance["sha256"][str(target.relative_to(snapshot))] = hashlib.sha256(target.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + name + "\n" for name, digest in sorted(provenance["sha256"].items())))
    print("Prepared equal-energy snapshot:", snapshot)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_equal_energy_target.py UNIQUE_SNAPSHOT_NAME")
    prepare(sys.argv[1])
