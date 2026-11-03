"""Plot four-model cloud energies and growth relative to independent linear growth.

Upper panels show absolute energies in Fa^2 M. Lower panels divide each model's
own energy growth E(t)/E(0) by exp(2 Im(omega) t/M), so initial domain-dependent
normalizations do not masquerade as changes in growth. Both this ratio and the
absolute E/E_linear ratio are exported. See export_growth_energies.jl for the
explicit, non-additive cloud-energy convention.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import plotstyle


def plot_energy_comparison(comparison, energy_prefix, output):
    hashes = {}

    def read(path, kind):
        path = Path(path)
        raw = path.read_bytes()
        hashes[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
        if kind == "csv":
            return np.atleast_1d(np.genfromtxt(io.BytesIO(raw), delimiter=",", names=True))
        if kind == "json":
            return json.loads(raw)
        return tomllib.loads(raw.decode())

    original = read(comparison, "json")
    post = read(str(energy_prefix) + "-postprocessing.toml", "toml")
    # Reject stale combinations of checkpoints, reports, trajectories and exports.
    expected_hashes = {**original["input_sha256"], **post["input_sha256"]}
    for path, expected in expected_hashes.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input has changed since the recorded comparison: {path}")
    for path in original["input_sha256"].keys() & post["input_sha256"].keys():
        if original["input_sha256"][path] != post["input_sha256"][path]:
            raise ValueError(f"Original comparison and postprocessing use different inputs: {path}")
    if not all(model["completed"] for model in original["models"].values()):
        raise ValueError("This plot requires completed input trajectories")
    rates = np.asarray(post["linear_amplitude_rates"])
    linear_e0 = np.asarray(post["linear_initial_energies"])
    if (rates.shape != (2,) or linear_e0.shape != (2,) or
            np.any(rates <= 0) or np.any(linear_e0 <= 0) or
            not np.all(np.isfinite(linear_e0)) or
            not np.allclose(rates, original["linear_rates"], rtol=1e-13, atol=0)):
        raise ValueError("Invalid linear reference")

    models = {}
    reports = {}
    for key in ("cosine", "chiral", "omiya"):
        if key == "omiya":
            paths = [p for p in original["input_sha256"] if p.endswith("-trajectory.csv")]
            if len(paths) != 1:
                raise ValueError("Require exactly one original Omiya trajectory")
            trajectory_path = paths[0]
            prefix = trajectory_path.removesuffix("-trajectory.csv")
        else:
            prefix = post[key]["prefix"]
        report = read(prefix + "-report.toml", "toml")
        reports[key] = report
        if (report["termination"] != "final_time" or
                report["alpha"] != original["alpha"] or report["spin"] != original["spin"] or
                report["potential"] != ("chiral" if key == "chiral" else "cosine")):
            raise ValueError(f"Incompatible {key} run")
        if key == "omiya":
            data = read(trajectory_path, "csv")
            energies = np.column_stack((data["primary_energy"], data["secondary_energy"]))
            steps = np.arange(len(data))
        else:
            path = str(energy_prefix) + f"-{key}.csv"
            data = read(path, "csv")
            if hashes[str(Path(path).resolve())] != post[key]["output_sha256"]:
                raise ValueError(f"Changed energy export: {key}")
            energies = np.column_stack((data["E211"], data["E322"]))
            steps = data["step"].astype(int)
            trajectory = read(prefix + f"-steps{report['steps']}.csv", "csv")
            if (not np.array_equal(data["time"], trajectory["time"][steps]) or
                    not np.array_equal(data["total_energy"], trajectory["energy"][steps]) or
                    not np.array_equal(data["A1"], trajectory["A1"][steps]) or
                    not np.array_equal(data["A2"], trajectory["A2"][steps])):
                raise ValueError(f"Energy export does not match saved trajectory: {key}")
        times = data["time"]
        if (len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0) or
                times[-1] != report["requested_final_time"] or
                not np.all(np.isfinite(energies)) or np.any(energies <= 0)):
            raise ValueError(f"Invalid energy trajectory: {key}")
        models[key] = {"times": times, "steps": steps, "energies": energies}
    times = models["omiya"]["times"]
    if any(m["times"][-1] != times[-1] for m in models.values()):
        raise ValueError("Mismatched end times")
    if reports["cosine"]["spatial_settings"] != reports["chiral"]["spatial_settings"]:
        raise ValueError("Joint domains differ")
    seeds = reports["cosine"]["initial_amplitudes"]
    if any(not np.allclose(r["initial_amplitudes"], seeds, rtol=1e-12, atol=0)
           for r in reports.values()):
        raise ValueError("Mismatched amplitude seeds")
    models["linear"] = {"times": times, "steps": np.arange(len(times)),
                        "energies": linear_e0 * np.exp(2 * times[:, None] * rates)}

    order = ("linear", "omiya", "cosine", "chiral")
    labels = {"linear": "Linear", "omiya": "Omiya cosine",
              "cosine": "Joint cosine", "chiral": "Joint chiral"}
    styles = {"linear": dict(color="#555555", linestyle=":", linewidth=2.0),
              "omiya": dict(color=plotstyle.PALETTE[1], linestyle="--", linewidth=2.3),
              "cosine": dict(color=plotstyle.PALETTE[2], linestyle="-", linewidth=1.8),
              "chiral": dict(color=plotstyle.PALETTE[0], linestyle="-.", linewidth=2.0)}
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), sharex=True,
                                gridspec_kw={"height_ratios": [1.3, 1]})
    figure.subplots_adjust(top=0.80, bottom=0.215, hspace=0.20, wspace=0.30)
    handles, summary, rows = {}, {}, []
    for key in order:
        model = models[key]
        t, energies = model["times"], model["energies"]
        exponential = np.exp(2 * t[:, None] * rates)
        growth = energies / energies[0]
        ratio = growth / exponential
        absolute_ratio = energies / (linear_e0 * exponential)
        for cloud in range(2):
            handles[key], = axes[0, cloud].plot(t/1e7, energies[:, cloud], **styles[key])
            axes[1, cloud].plot(t/1e7, ratio[:, cloud], **styles[key])
        summary[key] = {"samples": len(t), "completed": True, "last_time": float(t[-1]),
                        "initial_energies": energies[0].tolist(), "final_energies": energies[-1].tolist(),
                        "final_growth_factors": growth[-1].tolist(),
                        "final_growth_ratios_to_linear": ratio[-1].tolist(),
                        "final_absolute_ratios_to_linear": absolute_ratio[-1].tolist()}
        for index in range(len(t)):
            rows.append([key, model["steps"][index], t[index], *energies[index],
                         *growth[index], *ratio[index], *absolute_ratio[index]])
    for cloud, label in enumerate(("211", "322")):
        axes[0, cloud].set_title(label + " cloud", fontsize=13, pad=10)
        axes[0, cloud].set_ylabel(r"$E_{%s}/(F_a^2 M)$" % label, fontsize=12)
        axes[1, cloud].set_ylabel("Energy growth / linear growth", fontsize=10.5)
        axes[1, cloud].set_xlabel(r"$t/(10^7 M)$", fontsize=12)
        axes[1, cloud].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[0, 0].set_yscale("log")
    axes[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=False)
    for axis in axes.flat:
        axis.set_xlim(0, times[-1]/1e7)
        axis.grid(alpha=0.18)
        axis.tick_params(labelsize=10)
    figure.suptitle("Two-cloud energy evolution", fontsize=18, y=0.975)
    figure.text(0.5, 0.928,
                rf"Fixed Kerr background: $\alpha={original['alpha']:.2f}$, $\chi={original['spin']:.2f}$"
                + r"; $A_{211}(0)=10^{-6}$, $A_{322}(0)=10^{-7}$", ha="center", fontsize=11)
    figure.legend([handles[k] for k in order], [labels[k] for k in order],
                  loc="upper center", bbox_to_anchor=(0.5, 0.905), ncol=4, frameon=False, fontsize=11)
    joint = reports["cosine"]["spatial_settings"]
    omiya = reports["omiya"]["spatial_settings"]
    notes = [
        r"Lower panels: $[E_i(t)/E_i(0)]/\exp[2\,\mathrm{Im}(\omega_i)t]$; each model uses its own initial energy.",
        r"$E_{211}$: isolated nonlinear primary. $E_{322}$: quadratic secondary with averaged $U''(\Phi_{211})$.",
        "Mixed modes excluded from cloud labels; these energies are not a unique split of total energy.",
        f"Preliminary fixed-BH comparison · Joint/linear domain [{joint['inner']:g}, {joint['outer']:g}]; "
        f"Omiya [{omiya['inner']:g}, {omiya['outer']:g}] in r*/M.",
        "Matched outer amplitudes, not matched initial energies · All trajectories complete (1,000 steps).",
    ]
    for y, note in zip((0.145, 0.113, 0.087, 0.055, 0.029), notes):
        figure.text(0.5, y, note, ha="center", fontsize=8.8, color="#444444")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(str(output) + "." + extension, dpi=200)
    plt.close(figure)
    with Path(str(output) + "-trajectories.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["model", "step", "time", "E211", "E322", "growth211", "growth322",
                         "growth_ratio211", "growth_ratio322", "absolute_ratio211", "absolute_ratio322"])
        writer.writerows(rows)
    metadata = {"energy_convention": post["convention"], "energy_unit": post["energy_unit"],
                "lower_panels": "[E_i(t)/E_i(0)]/exp(2*Im(omega_i)*t/M)",
                "absolute_ratio_columns": "E_i(t)/[E_linear_i(0)*exp(2*Im(omega_i)*t/M)]",
                "linear_initial_energies": linear_e0.tolist(), "linear_amplitude_rates": rates.tolist(),
                "normalization": post["normalization"], "models": summary,
                "alpha": original["alpha"], "spin": original["spin"],
                "domains": {key: r["spatial_settings"] for key, r in reports.items()},
                "postprocessing": str(Path(str(energy_prefix) + "-postprocessing.toml").resolve()),
                "original_comparison": str(Path(comparison).resolve()), "input_sha256": hashes,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(str(output) + "-comparison.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path, help="Original amplitude comparison JSON")
    parser.add_argument("--energy-prefix", required=True, type=Path, help="Julia energy export prefix")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plot_energy_comparison(args.comparison, args.energy_prefix, args.output)
