"""TODO: describe verification of files recorded in a campaign source snapshot."""

import hashlib
import json
from pathlib import Path
import sys


def verify_snapshot(directory):
    """TODO: document the manifest format and failure on missing or changed files."""
    directory = Path(directory)
    manifest = json.loads((directory / "provenance.json").read_text())
    for name, expected in manifest["sha256"].items():
        path = directory / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError("source snapshot changed: " + name)
    print("Verified {} snapshot files".format(len(manifest["sha256"])))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_snapshot.py SNAPSHOT_DIRECTORY")
    verify_snapshot(sys.argv[1])
