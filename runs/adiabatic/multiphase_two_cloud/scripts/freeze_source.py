"""TODO: describe the campaign snapshot and its provenance files."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def freeze_source(name, seed_path=None):
    """TODO: document the snapshot name, destination, and no-overwrite behavior."""
    root = Path(__file__).resolve().parents[4]
    campaign = Path("runs/adiabatic/multiphase_two_cloud")
    if not name or Path(name).name != name or name in (".", ".."):
        raise ValueError("use a single directory name for the snapshot")
    destination = root / campaign / "work" / name
    destination.mkdir(parents=True, exist_ok=False)
    for relative in (Path("src"), Path("tests"), campaign / "scripts"):
        shutil.copytree(root / relative, destination / relative,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for name in ("Project.toml", "Manifest.toml"):
        shutil.copy2(root / name, destination / name)
    if seed_path is not None:
        (destination / "inputs").mkdir()
        shutil.copy2(seed_path, destination / "inputs" / "initial.jls")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, universal_newlines=True).strip()
    status = subprocess.check_output(["git", "status", "--short"], cwd=root, universal_newlines=True)
    diff = subprocess.check_output(["git", "diff", "HEAD", "--binary"], cwd=root)
    (destination / "working-tree.patch").write_bytes(diff)
    (destination / "working-tree-status.txt").write_text(status)
    hashes = {str(path.relative_to(destination)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(destination.rglob("*")) if path.is_file()}
    (destination / "provenance.json").write_text(json.dumps(
        {"git_revision": revision, "source_root": str(root), "sha256": hashes}, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: freeze_source.py SNAPSHOT_NAME [SEED_JLS]")
    freeze_source(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else None)
