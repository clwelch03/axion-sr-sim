"""Freeze the multimode solvers and seeds for energy-matched Omiya references."""
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
    parser.add_argument("result_prefixes", type=Path, nargs="+")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    campaign = root / "runs/adiabatic/multiphase_two_cloud"
    cases = []
    for prefix in args.result_prefixes:
        with Path(str(prefix) + ".toml").open("rb") as stream:
            metadata = tomllib.load(stream)
        if not metadata["passed"]:
            raise ValueError(f"Unaccepted result: {prefix}")
        cases.append({"name": prefix.name, **{key: metadata[key] for key in
            ("alpha", "spin", "inner", "outer", "target_energies")}})
    freeze_source(args.name)
    snapshot = campaign / "work" / args.name
    inputs = snapshot / "inputs"
    inputs.mkdir()
    seeds = campaign / "work/direct_fig7_20261026_v5/inputs"
    for source, destination in (("alpha030_cosine_branch.jls", "primary.jls"),
                                ("alpha030_cosine_secondary.jls", "secondary.jls")):
        shutil.copy2(seeds / source, inputs / destination)
    with (campaign / "results/equal_energy_20261028_v1.toml").open("rb") as stream:
        reference_energy = tomllib.load(stream)["reference_energy"]
    lines = [f"reference_energy = {reference_energy!r}"]
    for case, prefix in zip(cases, args.result_prefixes):
        lines.append("[[cases]]")
        lines.extend(f"{key} = {json.dumps(value)}" for key, value in case.items())
        shutil.copy2(Path(str(prefix) + ".toml"), inputs / (prefix.name + "-joint.toml"))
    (inputs / "cases.toml").write_text("\n".join(lines) + "\n")
    provenance = json.loads((snapshot / "provenance.json").read_text())
    for path in inputs.iterdir():
        provenance["sha256"][str(path.relative_to(snapshot))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (snapshot / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (snapshot / "source-manifest.sha256").write_text("".join(
        digest + "  " + name + "\n" for name, digest in sorted(provenance["sha256"].items())))


if __name__ == "__main__":
    main()
