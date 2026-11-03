"""Freeze the linear reference grid or the five missing secondary retries."""
import argparse
import hashlib
import json
from pathlib import Path
import tomllib
from freeze_source import freeze_source


def prepare(name, retry_missing_secondaries=False):
    root = Path(__file__).resolve().parents[4]
    freeze_source(name)
    snapshot = root / "runs/adiabatic/multiphase_two_cloud/work" / name
    (snapshot / "inputs").mkdir()
    lines = ['normalization = "first outer coefficient = 1"', ""]
    if retry_missing_secondaries:
        lines.extend(['clouds = ["secondary"]', 'frequency_precision_tolerance = true',
                      'coarse_tolerance = 1e-5', ""])
    for alpha in (0.3, 0.2, 0.4):
        for spin in (0.99, 0.95, 0.9, 0.8):
            if retry_missing_secondaries and not (alpha == 0.2 or (alpha == 0.3 and spin == 0.9)):
                continue
            case_name = f"alpha{round(100*alpha):03d}_spin{round(100*spin):03d}"
            lines.extend(["[[cases]]", f'name = "{case_name}"', f"alpha = {alpha!r}",
                          f"spin = {spin!r}", f"outer = {675.0 if alpha == 0.2 else 300.0}", ""])
    config_path = snapshot / "inputs/linear_grid.toml"
    config_path.write_text("\n".join(lines))
    count = 5 if retry_missing_secondaries else 12
    assert len(tomllib.loads(config_path.read_text())["cases"]) == count
    provenance = json.loads((snapshot / "provenance.json").read_text())
    provenance["sha256"]["inputs/linear_grid.toml"] = hashlib.sha256(config_path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + path + "\n" for path, digest in sorted(provenance["sha256"].items())))
    print(f"Prepared {count} linear reference backgrounds:", config_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name", help="fresh snapshot directory name")
    parser.add_argument("--retry-missing-secondaries", action="store_true",
                        help="retry alpha=0.2 at all four spins and alpha=0.3 at spin=0.9")
    args = parser.parse_args()
    prepare(args.name, args.retry_missing_secondaries)
