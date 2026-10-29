"""Plot six physical profiles and signed deformations at fixed outer amplitudes.

Usage: PYTHONPATH=src python plot_collocation_comparison.py RESULT_PREFIX FIGURE_PREFIX
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import plotstyle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_prefix", type=Path)
    parser.add_argument("figure_prefix", type=Path)
    parser.add_argument("--energy-references", type=Path,
                        help="TOML metadata for the equal-energy linear/Omiya reference export")
    args = parser.parse_args()
    metadata_path = Path(str(args.result_prefix) + ".toml")
    data_path = Path(str(args.result_prefix) + ".csv")
    with metadata_path.open("rb") as stream:
        metadata = tomllib.load(stream)
    equal_energies = metadata["normalization"] == "equal_cloud_energies"
    if not metadata["passed"] or metadata["normalization"] not in (
            "fixed_outer_coefficients", "equal_cloud_energies"):
        raise ValueError("Require an assessed result with a known normalization")
    if equal_energies and args.energy_references is not None:
        raise ValueError("Equal-energy exports already contain normalized references")
    if equal_energies:
        measured = np.asarray(metadata["joint_cloud_energies"])
        if np.max(np.abs(measured / metadata["reference_energy"] - 1)) >= 1e-4:
            raise ValueError("Joint cloud energies miss the requested comparison target")
    names = [model + "_" + cloud for model in ("linear", "omiya", "joint")
             for cloud in ("primary", "secondary")]
    curves = {name: [] for name in names}
    with data_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            curves[row["curve"]].append((float(row["r_star_over_M"]), float(row["R_abs"])))
    arrays = {name: np.asarray(values) for name, values in curves.items()}
    positions = arrays["linear_primary"][:, 0]
    for name, values in arrays.items():
        if values.shape != (len(positions), 2) or not np.isfinite(values).all():
            raise ValueError("Invalid curve: " + name)
        if not np.array_equal(values[:, 0], positions) or np.any(np.diff(positions) <= 0):
            raise ValueError("Require one shared increasing radial grid")
    energy_metadata = None
    amplitude_factors = {name: 1.0 for name in names}
    if args.energy_references is not None:
        with args.energy_references.open("rb") as stream:
            energy_metadata = tomllib.load(stream)
        if energy_metadata["normalization"] != "equal_cloud_energies":
            raise ValueError("Expected the verified equal-energy reference export")
        if energy_metadata["omiya_outer_amplitudes"] != metadata["cloud_amplitudes"]:
            raise ValueError("Omiya amplitudes differ from the joint run")
        if any(energy_metadata[key] != metadata[key] for key in ("spin", "alpha", "inner", "outer")):
            raise ValueError("Reference and joint settings differ")
        original_hash = metadata["reference_curve_metadata"]["input_sha256"]["direct_fig7_20261026_v5-input.jls"]
        if energy_metadata["input_sha256"]["direct_fig7_20261026_v5-input.jls"] != original_hash:
            raise ValueError("Reference exports use different joint inputs")
        for name, factor in energy_metadata["linear_scale_factors"].items():
            amplitude_factors[name] = factor
            arrays[name][:, 1] *= factor
    for suffix in (".png", ".pdf", "-plot-data.csv", "-curves.csv", "-plot.json"):
        if Path(str(args.figure_prefix) + suffix).exists():
            raise FileExistsError(str(args.figure_prefix) + suffix)
    args.figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(10.5, 8.0), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=(2, 1.2))
    upper = figure.add_subplot(grid[0, :])
    lower = [figure.add_subplot(grid[1, index]) for index in range(2)]
    plotted = {"r_star_over_M": positions}
    peaks = {}
    domain = (positions >= -100) & (positions <= 200)
    styles = {"linear": ":", "omiya": "--", "joint": "-"}
    for index, cloud in enumerate(("primary", "secondary")):
        color = plotstyle.PALETTE[1 if index == 0 else 4]
        linear = arrays["linear_" + cloud][:, 1]
        peak = float(np.max(linear[domain]))
        if peak <= 0:
            raise ValueError("Require a nonzero linear reference peak")
        peaks[cloud] = peak
        differences = {}
        for model in ("linear", "omiya", "joint"):
            name = model + "_" + cloud
            values = arrays[name][:, 1]
            plotted[name] = values
            upper.plot(positions, values, color=color, linestyle=styles[model],
                       linewidth=2, label=f"{cloud.capitalize()}: {model.capitalize()}")
            differences[model] = (values - linear) / peak
            if model != "linear":
                plotted[cloud + "_D_" + model] = differences[model]
                lower[index].plot(positions, 100 * differences[model], color=color,
                                  linestyle=styles[model], linewidth=1.8, label=model.capitalize())
        extra = differences["joint"] - differences["omiya"]
        plotted[cloud + "_extra"] = extra
        lower[index].plot(positions, 100 * extra, color="0.2", linewidth=1.3,
                          label="Joint − Omiya")
        lower[index].axhline(0, color="0.6", linestyle=":", linewidth=0.8)
        lower[index].set_title(cloud.capitalize())
        lower[index].set_ylabel("Signed change / linear peak [%]")
        lower[index].set_xlabel(r"$r_*/M$")
        lower[index].legend(frameon=False, fontsize=8)
    upper.set_ylabel(r"$|R_{lm\omega}|$")
    upper.set_xlabel(r"$r_*/M$")
    upper.set_ylim(bottom=0)
    upper.legend(frameon=False, ncol=2, fontsize=9)
    for axis in (upper, *lower):
        axis.set_xlim(-100, 200)
        axis.grid(alpha=0.18)
        axis.tick_params(direction="in", top=True, right=True)
    first_amplitude, second_amplitude = metadata["cloud_amplitudes"]
    normalization_caption = (
        f"Linear/Omiya reference E={energy_metadata['reference_energy']:.6g}; joint uses Omiya outer amplitudes"
        if energy_metadata is not None else
        f"Fixed outer amplitudes: A₁={first_amplitude:.8g}, A₂={second_amplitude:.8g}")
    if equal_energies:
        normalization_caption = (
            f"All clouds E₁=E₂≈{metadata['reference_energy']:.6g} Fₐ²M (Omiya energy convention)")
    figure.suptitle(
        rf"Preliminary joint-deformation comparison: $\alpha={metadata['alpha']:g}$, "
        rf"$a/M={metadata['spin']:g}$" + "\n"
        + normalization_caption + "\n"
        + f"Joint defect={metadata['independent_off_mesh_defect']:.2g}; {metadata['mesh_points']} mesh points",
        fontsize=11)
    for suffix in (".png", ".pdf"):
        output = Path(str(args.figure_prefix) + suffix)
        figure.savefig(output, dpi=200, bbox_inches="tight")
        print(output)
    plt.close(figure)
    with Path(str(args.figure_prefix) + "-plot-data.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(plotted)
        writer.writerows(zip(*plotted.values()))
    with data_path.open(newline="") as source, Path(str(args.figure_prefix) + "-curves.csv").open("w", newline="") as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            factor = amplitude_factors[row["curve"]]
            for key in ("R_real", "R_imag", "R_abs", "Q_real", "Q_imag"):
                row[key] = format(float(row[key]) * factor, ".17g")
            writer.writerow(row)
    hashed_paths = [metadata_path, data_path, Path(__file__)]
    if args.energy_references is not None:
        hashed_paths.append(args.energy_references)
    report = {"normalization": "equal_cloud_energies" if equal_energies else (
              "energy_normalized_references_joint_fixed_outer_coefficients"
              if energy_metadata is not None else "fixed_outer_coefficients"), "linear_peaks": peaks,
              "reference_energy_metadata": energy_metadata, "amplitude_factors": amplitude_factors,
              "deformation": "D_X=(abs(R_X)-abs(R_linear))/linear_peak; extra=D_joint-D_omiya",
              "plot_domain": [-100, 200], "saved_data_domain": [float(positions[0]), float(positions[-1])],
              "result_metadata": metadata, "input_sha256": {
                  str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in hashed_paths},
              "validation_scope": "Finite-domain radial check; angular, boundary, and retained-support systematics remain."}
    Path(str(args.figure_prefix) + "-plot.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
