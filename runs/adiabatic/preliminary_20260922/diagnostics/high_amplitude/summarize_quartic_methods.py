"""Summarize isolated method runs; never modify production checkpoints."""
import csv
import math
from pathlib import Path

directory = Path(__file__).resolve().parent / "quartic_methods"
directory.mkdir(exist_ok=True)
rows = []
for path in sorted(directory.glob("*.csv")):
    with path.open() as stream:
        rows.extend(csv.DictReader(stream))
rows.sort(key=lambda row: float(row["gain_per_minute"]), reverse=True)
lines = ["# Quartic continuation method comparison", "",
         "All cases start at A0=5.5837601798 and target 5.7 with the same grids,",
         "1e-9 residual tolerance, 12-correction attempt limit, and 600-second",
         "budget checked between attempts. SLURM allows 15 minutes per case.", "",
         "| Method | Last A0 | Gain/min | Seconds | Accepted | Rejected | Status |",
         "|---|---:|---:|---:|---:|---:|---|"]
for row in rows:
    lines.append(f"| {row['label']} | {float(row['amplitude']):.8f} | "
                 f"{float(row['gain_per_minute']):.5g} | {float(row['elapsed']):.1f} | "
                 f"{row['accepted']} | {row['rejected']} | {row['status']} |")
lines.extend(["", f"Results available for {len(rows)} of 7 cases.", ""])
completed = [row for row in rows if row["status"] == "target_reached"]
if len(completed) >= 2:
    reference = completed[0]
    tolerances = {"energy": 1e-7, "omega_real": 1e-9, "omega_imag": 1e-5, "total_flux": 1e-5}
    agreement = all(math.isclose(float(row[key]), float(reference[key]), rel_tol=tolerance)
                    for row in completed[1:] for key, tolerance in tolerances.items())
    lines.append(f"Common-endpoint energy, frequency, and flux agreement: {'PASS' if agreement else 'FAIL'}.")
    if agreement:
        lines.append(f"Fastest completed case: **{reference['label']}**. This is the candidate for a further production trial.")
    else:
        lines.append("Investigate endpoint disagreement before selecting a method.")
elif rows:
    lines.append(f"Highest measured progress rate: **{rows[0]['label']}**. Fewer than two cases reached the common endpoint; the ranking is provisional.")
else:
    lines.append("No result files were produced; inspect the array logs.")
lines.extend(["", "No production settings or checkpoints were changed by these experiments.", ""])
report = directory / "SUMMARY.md"
report.write_text("\n".join(lines))
print(report)
print("\n".join(lines))
