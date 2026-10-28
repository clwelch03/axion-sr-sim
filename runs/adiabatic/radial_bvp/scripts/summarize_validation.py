"""Extract numerical validation outputs from the preserved Julia job log."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

log_path, output_path = map(Path, sys.argv[1:])
rows = []
for line in log_path.read_text().splitlines():
    if line.startswith("FILE "):
        rows.append({"file": line[5:]})
    elif line.startswith("finite_domain_observables_h"):
        label, content = line.split("=", 1)
        rows[-1]["observable_grid_spacing"] = float(label.split("_h")[-1])
        rows[-1]["observables"] = {
            key: float(value)
            for key, value in re.findall(r"(\w+) = ([-+0-9.eE]+)", content)
        }
    elif line.startswith("quadrature_relative_changes_h1_to_h0.5="):
        rows[-1]["quadrature_relative_changes"] = ast.literal_eval(line.split("=", 1)[1])
    elif line.startswith("short_reintegration="):
        for key, values in re.findall(r"(defects|step_changes) = (\[[^\]]*\])", line):
            rows[-1]["short_reintegration_" + key] = ast.literal_eval(values)
for row in rows:
    assert len(row["observables"]) == 6
    assert len(row["short_reintegration_defects"]) == 28
    assert len(row["short_reintegration_step_changes"]) == 28
assert rows, "No completed validation records found"
output_path.write_text(json.dumps({
    "log": str(log_path), "log_sha256": hashlib.sha256(log_path.read_bytes()).hexdigest(),
    "runs": rows,
}, indent=2) + "\n")
