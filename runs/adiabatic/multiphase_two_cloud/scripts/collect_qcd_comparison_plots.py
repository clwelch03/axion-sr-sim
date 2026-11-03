"""Collect completed QCD comparisons into a 36-point manifest and PDF bundle.

Do not label the bundle complete unless every expected point has passed and
has its own three-model figure. Missing cases remain explicit in the manifest.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tomllib
from pypdf import PdfWriter

ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / "src/common/plotstyle.py").is_file())


def collect(campaign, output):
    output.mkdir(parents=True, exist_ok=True)
    with (campaign / "results/qcd_populations_20261029_v1/population_targets.csv").open(newline="") as stream:
        populations = list(csv.DictReader(stream))
    rows, files, writer = [], {}, PdfWriter()
    for population in populations:
        alpha, spin, mass = (float(population[key]) for key in ("alpha", "spin", "mass_solar"))
        mode = population["dominant_mode"]
        name = f"qcd{mode}_alpha{round(100*alpha):03d}_spin{round(100*spin):03d}_mass{round(mass):03d}"
        if (alpha, spin, mass) == (.3, .99, 5.):
            directory = campaign / "results/qcd_joint_probe_20261030_v2"
        else:
            version = "qcd_211_comparisons_20261030_v2" if mode == "211" else "qcd_322_comparisons_20261030_v1"
            directory = campaign / "results" / version / name
        row = dict(alpha=alpha, spin=spin, mass_solar=mass, dominant_mode=mode,
                   status="pending_or_failed", figure=str(directory / "both_clouds.pdf"))
        report_path, plot_path = directory / "comparison.toml", directory / "plot_metadata.json"
        if report_path.exists() and plot_path.exists() and (directory / "both_clouds.pdf").exists():
            report = tomllib.loads(report_path.read_text())
            plot = json.loads(plot_path.read_text())
            assert report["passed"] and report["joint_computed"]
            assert set(plot["models"]) == {"linear", "omiya", "joint"}
            assert report["secondary_unit_matching_norm"] <= report["secondary_matching_tolerance"]
            assert (report["case"]["alpha"], report["case"]["spin"], report["case"]["mass_solar"]) == (alpha, spin, mass)
            row["status"] = "plotted"
            for cloud in ("211", "322"):
                row[f"energy_{cloud}"] = report[f"energy_{cloud}"]
                measurements = plot["clouds"][cloud]
                row[f"omiya_linear_peak_percent_{cloud}"] = measurements["maximum_difference_over_linear_peak_percent"]
                row[f"joint_omiya_peak_percent_{cloud}"] = measurements["maximum_joint_minus_omiya_over_linear_peak_percent"]
            writer.append(str(directory / "both_clouds.pdf"))
            for path in (report_path, plot_path, directory / "both_clouds.pdf"):
                files[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(row)
    assert len(rows) == 36
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (output / "comparison_manifest.csv").open("w", newline="") as stream:
        table = csv.DictWriter(stream, fieldnames=columns)
        table.writeheader()
        table.writerows(rows)
    completed = sum(row["status"] == "plotted" for row in rows)
    bundle = "all_36_profiles.pdf" if completed == 36 else f"completed_{completed}_of_36_profiles.pdf"
    if completed:
        writer.write(str(output / bundle))
    writer.close()
    (output / "metadata.json").write_text(json.dumps(dict(completed=completed, expected=36,
        complete=completed == 36, bundle=bundle if completed else None,
        potential="cosine", population_stage="original first-peak seed diagnostic",
        joint_regime="nonlinear dominant cloud plus coupled small-secondary limit",
        spatial_convergence_established=False, source_sha256=files), indent=2) + "\n")
    print(f"Collected {completed}/36 QCD comparisons into {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    collect(args.campaign.resolve(), args.output.resolve())
