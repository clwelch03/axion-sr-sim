"""Freeze existing linear solvers and the twelve requested alpha/spin pairs."""
import hashlib
import json
from pathlib import Path
import sys
import tomllib
from freeze_source import freeze_source


def prepare(name):
    root = Path(__file__).resolve().parents[4]
    freeze_source(name)
    snapshot = root / "runs/adiabatic/multiphase_two_cloud/work" / name
    (snapshot / "inputs").mkdir()
    lines = ['normalization = "first outer coefficient = 1"', ""]
    for alpha in (0.3, 0.2, 0.4):
        for spin in (0.99, 0.95, 0.9, 0.8):
            case_name = f"alpha{round(100*alpha):03d}_spin{round(100*spin):03d}"
            lines.extend(["[[cases]]", f'name = "{case_name}"', f"alpha = {alpha!r}",
                          f"spin = {spin!r}", f"outer = {675.0 if alpha == 0.2 else 300.0}", ""])
    config_path = snapshot / "inputs/linear_grid.toml"
    config_path.write_text("\n".join(lines))
    assert len(tomllib.loads(config_path.read_text())["cases"]) == 12
    provenance = json.loads((snapshot / "provenance.json").read_text())
    provenance["sha256"]["inputs/linear_grid.toml"] = hashlib.sha256(config_path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + path + "\n" for path, digest in sorted(provenance["sha256"].items())))
    print("Prepared 12 linear reference pairs:", config_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_linear_reference_grid.py UNIQUE_NAME")
    prepare(sys.argv[1])
