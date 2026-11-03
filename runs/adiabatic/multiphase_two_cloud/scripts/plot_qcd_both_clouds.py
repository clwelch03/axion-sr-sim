"""Plot both clouds in an energy-matched QCD comparison.

Each cloud has its own physical vertical scale; no curve is peak-normalized.
The difference panels use that cloud's linear peak only as a reporting unit.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / "src/common/plotstyle.py").is_file())
sys.path.insert(0, str(ROOT / "src"))
from common import plotstyle


def plot(prefix):
    report = tomllib.loads(prefix.with_suffix(".toml").read_text())
    assert report["passed"]
    dominant_mode = report["case"]["dominant_mode"]
    minor_mode = "322" if dominant_mode == "211" else "211"
    joint_computed = report["joint_computed"]
    models = ("linear", "omiya", "joint") if joint_computed else ("linear", "omiya")
    for key in (f"linear_energy_{minor_mode}", f"omiya_energy_{minor_mode}"):
        if key in report:
            assert abs(report[key] / report[f"energy_{minor_mode}"] - 1) < 1e-12
    with prefix.with_suffix(".csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    heights = [2, 1, 1] if joint_computed else [2, 1]
    fig, axes = plt.subplots(len(heights), 2, figsize=(11, 8.5 if joint_computed else 6.6), sharex="col",
        gridspec_kw={"height_ratios": heights, "hspace": .17, "wspace": .34})
    measurements = {}
    for column, (cloud, degree) in enumerate((("211", 1), ("322", 2))):
        curves = {}
        for model in models:
            selected = sorted((row for row in rows if row["cloud"] == cloud
                and row["model"] == model and int(row["ell"]) == int(row["m"]) == degree),
                key=lambda row: (float(row["r_star_over_M"]), row["branch"]))
            curves[model] = np.array([[float(row["r_star_over_M"]), float(row["R_abs"])]
                                     for row in selected])
            assert curves[model].size and np.all(np.isfinite(curves[model]))
        assert np.allclose(curves["linear"][:, 0], curves["omiya"][:, 0], rtol=0, atol=1e-9)
        positions, linear = curves["linear"].T
        omiya = curves["omiya"][:, 1]
        difference = 100 * (omiya - linear) / np.max(linear)
        upper, lower = axes[:2, column]
        if joint_computed:
            assert np.allclose(positions, curves["joint"][:, 0], rtol=0, atol=1e-9)
            joint = curves["joint"][:, 1]
            upper.plot(positions, joint, color=plotstyle.PALETTE[2], lw=2.8,
                       label="Joint (probe limit)")
        upper.plot(positions, omiya, color=plotstyle.PALETTE[1], lw=2, label="Omiya")
        upper.plot(positions, linear, color="0.15", ls="--", lw=1.5, label="Linear")
        upper.set_ylabel(rf"$|R_{{{degree}{degree}}}|$")
        energy = report[f"energy_{cloud}"]
        if 1e-3 <= energy < 1e4:
            energy_text = f"{energy:.5f}"
        else:
            mantissa, exponent = f"{energy:.4e}".split("e")
            energy_text = rf"{mantissa}\times10^{{{int(exponent)}}}"
        upper.set_title(f"{cloud}" + (" dominant cloud" if cloud == dominant_mode else " subdominant probe")
            + "\n" + rf"$E/(F_a^2 M)={energy_text}$", fontsize=11)
        upper.set_ylim(bottom=0)
        upper.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
        upper.legend(frameon=False, fontsize=9, loc="upper right")
        lower.plot(positions, difference, color=plotstyle.PALETTE[1], lw=1.8, label="Omiya")
        if joint_computed:
            joint_difference = 100 * (joint-linear)/np.max(linear)
            lower.plot(positions, joint_difference, color=plotstyle.PALETTE[2], lw=1.5,
                       ls="--", label="Joint")
            between_models = 100 * (joint-omiya)/np.max(linear)
            axes[2, column].plot(positions, between_models, color=plotstyle.PALETTE[2], lw=1.8)
            axes[2, column].axhline(0, color=".5", lw=.8)
            axes[2, column].set_ylabel("Joint − Omiya\n[% of linear peak]", fontsize=10)
            axes[2, column].ticklabel_format(axis="y", style="sci", scilimits=(-3,3), useMathText=True)
            if np.max(np.abs(between_models)) < 1e-10:
                axes[2, column].text(.98, .8, "Roundoff level", transform=axes[2,column].transAxes,
                                     ha="right", fontsize=10)
            elif cloud == minor_mode:
                recovery_floor = 100*report["omiya_profile_recovery_without_feedback"]
                if np.max(np.abs(between_models)) <= 10*recovery_floor:
                    axes[2,column].axhspan(-recovery_floor,recovery_floor,color=".7",alpha=.3)
                    axes[2,column].text(.98,.8,"Difference not resolved by\nconsistency check",
                        transform=axes[2,column].transAxes,ha="right",fontsize=9)
        lower.axhline(0, color=".5", lw=.8)
        lower.set_ylabel(("Model − linear" if joint_computed else "Omiya − linear")
                         + "\n[% of linear peak]", fontsize=10)
        axes[-1,column].set_xlabel(r"$r_*/M$")
        for axis in axes[:,column]:
            axis.set_xlim(-100, 300)
            axis.grid(alpha=.18)
            axis.tick_params(direction="in", top=True, right=True)
        measurements[cloud] = dict(energy=energy,
            maximum_difference_over_linear_peak_percent=float(np.max(np.abs(difference))),
            linear_peak=float(np.max(linear)), omiya_peak=float(np.max(omiya)))
        if joint_computed:
            measurements[cloud].update(
                maximum_joint_minus_linear_over_linear_peak_percent=float(np.max(np.abs(joint_difference))),
                maximum_joint_minus_omiya_over_linear_peak_percent=float(np.max(np.abs(between_models))),
                joint_peak=float(np.max(joint)))
            if cloud == minor_mode:
                measurements[cloud]["consistency_check_floor_percent"] = 100*report["omiya_profile_recovery_without_feedback"]
                measurements[cloud]["difference_exceeds_ten_times_consistency_floor"] = bool(
                    np.max(np.abs(between_models)) > 1000*report["omiya_profile_recovery_without_feedback"])
    case = report["case"]
    fig.suptitle(rf"$M={case['mass_solar']:g}M_\odot,\ \alpha={case['alpha']:g},\ a_*={case['spin']:g}$" +
        "\nSame energy across models for each cloud; separate physical amplitude scales", y=.99)
    fig.subplots_adjust(top=.86, bottom=.13)
    fig.text(.5, .035,
        ("Joint: nonlinear primary and coupled small-secondary response; preliminary spatial resolution.\n"
         if joint_computed else "Linear/Omiya preview; joint profiles pending.\n")
        + rf"{minor_mode} uses the retained seed estimate $N_{{{minor_mode}}}\simeq{report['population']['secondary_reference_occupation']:.3f}$;"
        + " its normalization is illustrative at microscopic occupation.",
        ha="center", fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(prefix.parent / f"both_clouds.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    sources = [prefix.with_suffix(suffix) for suffix in (".csv", ".toml")]
    archived_driver = prefix.parent / Path(__file__).name
    archived_driver.write_bytes(Path(__file__).read_bytes())
    sources.extend((archived_driver, ROOT / "src/common/plotstyle.py"))
    metadata = dict(clouds=measurements, models=list(models), joint_computed=joint_computed,
        amplitude_normalization="Physical amplitudes; separate axes, no peak rescaling.",
        source_sha256={(str(path.resolve().relative_to(ROOT)) if path.resolve().is_relative_to(ROOT)
                       else str(path.resolve())): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sources})
    (prefix.parent / "plot_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(measurements, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path)
    plot(parser.parse_args().prefix)
