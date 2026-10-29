"""Plot accepted sweep profiles against energy-matched linear/Omiya references.

Usage: PYTHONPATH=src python plot_parameter_energy_cases.py REFERENCE_PREFIX
       FIGURE_PREFIX RESULT_PREFIX [RESULT_PREFIX ...]
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


CLOUDS = ("primary", "secondary")
MODES = {(1, 1, 0): "primary", (2, 0, 1): "secondary"}
COLORS = (plotstyle.PALETTE[1], plotstyle.PALETTE[4])
DOMAIN = (-100, 200)


def read_metadata(prefix):
    with Path(str(prefix) + ".toml").open("rb") as stream:
        return tomllib.load(stream)


def read_profiles(prefix, references=False, model="linear"):
    curves = {cloud: [] for cloud in CLOUDS}
    with Path(str(prefix) + ".csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if references:
                cloud = {model + "_" + name: name for name in CLOUDS}.get(row["curve"])
            else:
                cloud = MODES.get(tuple(int(row[key]) for key in ("ell", "q1", "q2")))
            if cloud is not None:
                radius, real, imag, magnitude = (float(row[key]) for key in
                    ("r_star_over_M", "R_real", "R_imag", "R_abs"))
                if not np.isclose(abs(complex(real, imag)), magnitude, rtol=1e-12, atol=1e-15):
                    raise ValueError(f"Inconsistent complex profile in {prefix}")
                curves[cloud].append((radius, magnitude))
    curves = {cloud: np.asarray(values) for cloud, values in curves.items()}
    for cloud, values in curves.items():
        if (values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all()
                or np.any(np.diff(values[:, 0]) <= 0)):
            raise ValueError(f"Invalid samples: {prefix}, {cloud}")
    if not np.array_equal(curves["primary"][:, 0], curves["secondary"][:, 0]):
        raise ValueError(f"Cloud grids differ in {prefix}")
    return curves


def case_title(case):
    first, second = case["energy_factors"]
    return rf"$E_1/E_0={first:g},\quad E_2/E_0={second:g}$"


def draw_case(figure, slot, case, limits):
    grid = slot.subgridspec(2, 1, height_ratios=(2.1, 1), hspace=0.06)
    upper = figure.add_subplot(grid[0])
    lower = figure.add_subplot(grid[1], sharex=upper)
    for cloud, color in zip(CLOUDS, COLORS):
        data = case["curves"][cloud]
        upper.plot(case["positions"], data["joint"], color=color, linewidth=2,
                   label=f"{cloud.capitalize()}: joint")
        upper.plot(case["positions"], data["linear"], color=color, linewidth=1.8,
                   linestyle=":", label=f"{cloud.capitalize()}: linear")
        if "omiya" in data:
            upper.plot(case["positions"], data["omiya"], color=color, linewidth=1.8,
                       linestyle="--", label=f"{cloud.capitalize()}: Omiya")
            lower.plot(case["positions"], data["joint_minus_omiya_percent"], color=color, linewidth=1.8)
        else:
            lower.plot(case["positions"], data["deformation_percent"], color=color, linewidth=1.8)
    upper.set_title(case_title(case), fontsize=12, pad=9)
    upper.set_ylabel(r"$|R_{\ell m\omega}|$")
    upper.set_ylim(0, limits["amplitude"])
    upper.tick_params(labelbottom=False)
    lower.axhline(0, color="0.5", linestyle=":", linewidth=0.8)
    lower.set_ylim(*limits["deformation"])
    lower.set_ylabel(("Joint − Omiya /" if case["omiya_metadata"] is not None else "Change /")
                     + "\nlinear peak [%]")
    lower.set_xlabel(r"$r_*/M$")
    for axis in (upper, lower):
        axis.set_xlim(*DOMAIN)
        axis.grid(alpha=0.18)
        axis.tick_params(direction="in", top=True, right=True)
    return upper.get_legend_handles_labels()


def save_figure(figure, prefix):
    for suffix in (".png", ".pdf"):
        output = Path(str(prefix) + suffix)
        figure.savefig(output, dpi=180, bbox_inches="tight")
        print(output)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_prefix", type=Path,
                        help="Equal-energy export with linear profiles on a common radial grid")
    parser.add_argument("figure_prefix", type=Path)
    parser.add_argument("result_prefixes", type=Path, nargs="+")
    parser.add_argument("--omiya-directory", type=Path,
                        help="Accepted Omiya exports with matching case names and energy targets")
    args = parser.parse_args()
    reference_metadata = read_metadata(args.reference_prefix)
    if not reference_metadata["passed"] or reference_metadata["normalization"] != "equal_cloud_energies":
        raise ValueError("Require assessed, energy-normalized references")
    references = read_profiles(args.reference_prefix, references=True)
    reference_energies = reference_metadata["reference_curve_metadata"]["energies_by_curve"]
    reference_energy = reference_metadata["reference_energy"]
    cases = []
    for prefix in args.result_prefixes:
        metadata = read_metadata(prefix)
        if not metadata["passed"]:
            raise ValueError(f"Result not accepted: {prefix}")
        for key in ("alpha", "spin", "inner", "outer"):
            if metadata[key] != reference_metadata[key]:
                raise ValueError(f"Reference {key} differs from {prefix}")
        targets = np.asarray(metadata["target_energies"])
        measured = np.asarray(metadata["joint_cloud_energies"])
        if (targets.shape != (2,) or measured.shape != (2,) or
                not np.isfinite(targets).all() or not np.isfinite(measured).all() or
                np.any(targets <= 0) or np.max(np.abs(measured / targets - 1)) >= 1e-4):
            raise ValueError(f"Invalid energy targets or unaccepted energy errors: {prefix}")
        joint = read_profiles(prefix)
        positions = joint["primary"][:, 0]
        displayed = (positions >= DOMAIN[0]) & (positions <= DOMAIN[1])
        omiya_metadata = None
        if args.omiya_directory is not None:
            omiya_prefix = args.omiya_directory / prefix.name
            omiya_metadata = read_metadata(omiya_prefix)
            if not omiya_metadata["passed"] or omiya_metadata["model"] != "omiya_multimode_primary_secondary":
                raise ValueError(f"Require an accepted Omiya result: {omiya_prefix}")
            for key in ("alpha", "spin", "inner", "outer", "target_energies"):
                if omiya_metadata[key] != metadata[key]:
                    raise ValueError(f"Omiya {key} differs from {prefix}")
            omiya_energies = np.asarray(omiya_metadata["measured_energies"])
            if (not np.isfinite(omiya_energies).all() or
                    np.max(np.abs(omiya_energies / targets - 1)) >= 1e-4):
                raise ValueError(f"Omiya energy targets not met: {omiya_prefix}")
            omiya = read_profiles(omiya_prefix, references=True, model="omiya")
        curves, factors = {}, {}
        for index, cloud in enumerate(CLOUDS):
            if not np.array_equal(positions, references[cloud][:, 0]):
                raise ValueError("Require shared reference/result grid; no interpolation is applied")
            factor = float(np.sqrt(targets[index] / reference_energies["linear_" + cloud]))
            linear = references[cloud][:, 1] * factor
            peak = float(np.max(linear[displayed]))
            if peak <= 0:
                raise ValueError("Require positive linear reference peak")
            curves[cloud] = {"joint": joint[cloud][:, 1], "linear": linear,
                             "deformation_percent": 100 * (joint[cloud][:, 1] - linear) / peak}
            if omiya_metadata is not None:
                if not np.array_equal(positions, omiya[cloud][:, 0]):
                    raise ValueError("Require shared Omiya/result grid")
                curves[cloud]["omiya"] = omiya[cloud][:, 1]
                curves[cloud]["omiya_deformation_percent"] = 100 * (omiya[cloud][:, 1] - linear) / peak
                curves[cloud]["joint_minus_omiya_percent"] = 100 * (joint[cloud][:, 1] - omiya[cloud][:, 1]) / peak
            factors[cloud] = factor
        cases.append({"name": prefix.name, "metadata": metadata, "positions": positions,
                      "energy_factors": (targets / reference_energy).tolist(),
                      "linear_amplitude_factors": factors, "curves": curves,
                      "omiya_metadata": omiya_metadata})

    prefixes = [args.figure_prefix] + [Path(str(args.figure_prefix) + "-" + case["name"]) for case in cases]
    if len(set(prefixes)) != len(prefixes):
        raise ValueError("Duplicate output names")
    outputs = [Path(str(prefix) + suffix) for prefix in prefixes for suffix in (".png", ".pdf")]
    outputs += [Path(str(args.figure_prefix) + suffix) for suffix in ("-plot-data.csv", "-plot.json")]
    for output in outputs:
        if output.exists():
            raise FileExistsError(output)
    args.figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_amplitudes, all_changes = [], []
    for case in cases:
        mask = (case["positions"] >= DOMAIN[0]) & (case["positions"] <= DOMAIN[1])
        for curve in case["curves"].values():
            all_amplitudes.extend([curve["joint"][mask].max(), curve["linear"][mask].max()])
            if "omiya" in curve:
                all_amplitudes.append(curve["omiya"][mask].max())
            change = curve["joint_minus_omiya_percent"] if "omiya" in curve else curve["deformation_percent"]
            all_changes.extend([change[mask].min(), change[mask].max()])
    bottom, top = min(0, min(all_changes)), max(0, max(all_changes))
    pad = max(1, (top - bottom) * 0.09)
    limits = {"amplitude": 1.08 * max(all_amplitudes), "deformation": (bottom - pad, top + pad)}
    rows, columns = (len(cases) + 1) // 2, min(2, len(cases))
    figure = plt.figure(figsize=(9, 7) if len(cases) == 1 else
                         (6.4 * columns, 4.9 * rows + 1.1), layout="constrained")
    grid = figure.add_gridspec(rows, columns)
    for index, case in enumerate(cases):
        handles, labels = draw_case(figure, grid[index // columns, index % columns], case, limits)
    comparison_caption = ("linear and Omiya references matched to each energy target"
                          if args.omiya_directory is not None else
                          "linear references matched to each energy target")
    figure.suptitle(
        rf"Joint cloud profiles: $\alpha={reference_metadata['alpha']:g},\ a/M={reference_metadata['spin']:g}$"
        + "\n" + rf"$E_0={reference_energy:.3f}\,F_a^2 M$; " + comparison_caption,
        fontsize=14)
    if args.omiya_directory is not None:
        order = (0, 3, 2, 5, 1, 4)
        handles, labels = [handles[i] for i in order], [labels[i] for i in order]
    figure.legend(handles, labels, loc="outside lower center",
                  ncol=3 if args.omiya_directory is not None else 4, frameon=False, fontsize=10)
    save_figure(figure, args.figure_prefix)
    for case, prefix in zip(cases, prefixes[1:]):
        figure = plt.figure(figsize=(9, 7), layout="constrained")
        handles, labels = draw_case(figure, figure.add_gridspec(1, 1)[0], case, limits)
        figure.suptitle(
            rf"$\alpha={reference_metadata['alpha']:g},\ a/M={reference_metadata['spin']:g},\ E_0={reference_energy:.3f}\,F_a^2 M$"
            + "\n" + comparison_caption.capitalize(), fontsize=12)
        figure.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
        save_figure(figure, prefix)

    with outputs[-2].open("w", newline="") as stream:
        writer = csv.writer(stream)
        columns = ["joint", "linear", "deformation_percent"]
        if args.omiya_directory is not None:
            columns += ["omiya", "omiya_deformation_percent", "joint_minus_omiya_percent"]
        writer.writerow(["case", "cloud", "r_star_over_M", *columns])
        for case in cases:
            for cloud, curve in case["curves"].items():
                writer.writerows((case["name"], cloud, *values) for values in zip(
                    case["positions"], *(curve[column] for column in columns)))
    inputs = [Path(__file__)] + [Path(str(prefix) + suffix) for prefix in
        [args.reference_prefix, *args.result_prefixes] for suffix in (".csv", ".toml")]
    if args.omiya_directory is not None:
        inputs += [Path(str(args.omiya_directory / case["name"]) + suffix)
                   for case in cases for suffix in (".csv", ".toml")]
    report = {"cases": [{key: value for key, value in case.items() if key not in ("positions", "curves")}
                        for case in cases],
              "reference_prefix": str(args.reference_prefix), "reference_energy": reference_energy,
              "omiya_directory": str(args.omiya_directory) if args.omiya_directory is not None else None,
              "plotted_modes": {cloud: list(mode) for mode, cloud in MODES.items()},
              "plot_domain": list(DOMAIN), "linear_scaling": "sqrt(target_energy / saved_linear_energy)",
              "joint_scaling": "none", "deformation": "100 * (abs(R_joint) - abs(R_linear)) / max_displayed(abs(R_linear))",
              "lower_panels": "100 * (abs(R_joint) - abs(R_omiya)) / max_displayed(abs(R_linear))"
                  if args.omiya_directory is not None else "deformation",
              "validation_scope": "Accepted finite-domain profiles; spatial/angular/boundary convergence remains separate.",
              "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}
    outputs[-1].write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
