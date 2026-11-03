"""Export independent linear growth and compare Omiya, joint cosine and chiral.

Linear amplitudes obey dA_i/d(t/M) = Im(omega_i) A_i at fixed BH parameters.
The plot compares outer-amplitude growth, not matched physical occupations.
"""
import argparse
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


def plot_growth_comparison(cosine, chiral, omiya, output, allow_partial=False):
    inputs = {}

    def read(path, table=False):
        path = Path(path)
        data = path.read_bytes()
        inputs[str(path.resolve())] = hashlib.sha256(data).hexdigest()
        return (np.atleast_1d(np.genfromtxt(io.BytesIO(data), delimiter=",", names=True))
                if table else tomllib.loads(data.decode()))

    prefixes = {"cosine": str(cosine), "chiral": str(chiral), "omiya": str(omiya)}
    reports = {key: read(prefix + "-report.toml") for key, prefix in prefixes.items()}
    reference = reports["cosine"]
    times = np.linspace(0, reference["requested_final_time"], reference["steps"] + 1)
    seeds = np.asarray(reference["initial_amplitudes"])
    models = {}
    for key, prefix in prefixes.items():
        report = reports[key]
        if report["potential"] != ("chiral" if key == "chiral" else "cosine"):
            raise ValueError("Wrong potential for " + key)
        for field in ("alpha", "spin", "requested_final_time"):
            if report[field] != reference[field]:
                raise ValueError("Mismatched " + field + " for " + key)
        if not np.allclose(report["initial_amplitudes"], seeds, rtol=1e-12, atol=0):
            raise ValueError("Mismatched seeds for " + key)
        if key != "omiya" and report["spatial_settings"] != reference["spatial_settings"]:
            raise ValueError("Joint spatial settings differ")
        suffix = "-trajectory.csv" if key == "omiya" else "-steps%d.csv" % reference["steps"]
        data = read(prefix + suffix, table=True)
        amplitudes = np.column_stack((data["A1"], data["A2"]))
        if (len(data) > len(times) or len(data) < 2 or
                not np.allclose(data["time"], times[:len(data)], rtol=1e-13, atol=1e-8) or
                not np.all(np.isfinite(amplitudes)) or np.any(amplitudes <= 0) or
                not np.allclose(amplitudes[0], seeds, rtol=1e-12, atol=0)):
            raise ValueError("Invalid trajectory or incompatible time grid for " + key)
        complete = len(data) == len(times) and report["termination"] == "final_time"
        if not complete and not allow_partial:
            raise ValueError(key + " has not completed; use --allow-partial for a labeled preview")
        if report["termination"] == "failed":
            raise ValueError(key + " reports a failed evolution")
        models[key] = {"time": data["time"], "amplitudes": amplitudes,
                       "complete": complete, "steps": len(data)-1}

    if models["cosine"]["complete"]:
        times = models["cosine"]["time"].copy()
    linear_reference = read(str(cosine) + "-linear-reference.csv", table=True)
    if len(linear_reference) != 2 or not np.array_equal(linear_reference["cloud"], [1, 2]):
        raise ValueError("Require linear modes in cloud order 211, 322")
    rates = linear_reference["linear_growth_rate"]
    if (not np.all(np.isfinite(rates)) or np.any(rates <= 0) or
            not np.allclose(linear_reference["initial_amplitude"], seeds, rtol=1e-12, atol=0) or
            not np.allclose(rates, reference["linear_amplitude_rates"], rtol=1e-13, atol=0) or
            not np.allclose(rates, reports["chiral"]["linear_amplitude_rates"], rtol=1e-13, atol=0)):
        raise ValueError("Inconsistent linear rates or normalization")
    linear = seeds * np.exp(times[:, None] * rates)
    models["linear"] = {"time": times, "amplitudes": linear,
                        "complete": True, "steps": len(times)-1}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(str(output) + "-linear.csv",
               np.column_stack((times, linear, np.broadcast_to(rates, linear.shape))),
               delimiter=",", header="time,A1,A2,log_rate1,log_rate2", comments="")

    labels = {"linear": "Linear (no interactions)", "omiya": "Omiya (cosine)",
              "cosine": "Joint cosine", "chiral": "Joint chiral"}
    styles = {"linear": dict(color="#555555", linestyle=":", linewidth=1.7),
              "omiya": dict(color=plotstyle.PALETTE[1], linestyle="--", linewidth=2.3),
              "cosine": dict(color=plotstyle.PALETTE[2], linestyle="-", linewidth=1.8),
              "chiral": dict(color=plotstyle.PALETTE[0], linestyle="-", linewidth=2.0)}
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.8), sharex=True)
    figure.subplots_adjust(top=0.80, bottom=0.16, hspace=0.28, wspace=0.27)
    handles = {}
    columns, names = [times], ["time"]
    summary = {}
    for key in ("linear", "cosine", "omiya", "chiral"):
        model = models[key]
        size = len(model["time"])
        growth = model["amplitudes"] / seeds
        deviations = 100 * (model["amplitudes"] / linear[:size] - 1)
        for cloud in range(2):
            handles[key], = axes[0, cloud].plot(model["time"]/1e7, growth[:, cloud],
                                               label=labels[key], **styles[key])
            axes[1, cloud].plot(model["time"]/1e7, deviations[:, cloud], **styles[key])
            padded = np.full(len(times), np.nan)
            padded[:size] = model["amplitudes"][:, cloud]
            columns.append(padded)
            names.append("%s_A%d" % (key, cloud+1))
        if not model["complete"]:
            for cloud in range(2):
                axes[0, cloud].plot(model["time"][-1]/1e7, growth[-1, cloud], "o",
                                    color=styles[key]["color"], markersize=4)
        summary[key] = {"completed": model["complete"], "steps": model["steps"],
                        "last_time": float(model["time"][-1]),
                        "last_amplitudes": model["amplitudes"][-1].tolist(),
                        "amplitude_growth_factors": growth[-1].tolist(),
                        "deviation_from_linear_percent": deviations[-1].tolist()}
    for cloud, label in enumerate(("211", "322")):
        axes[0, cloud].set_title(label + " cloud", fontsize=13)
        axes[0, cloud].set_ylabel(r"$A_{%s}(t)/A_{%s}(0)$" % (label, label))
        axes[1, cloud].set_ylabel("Deviation from linear (%)")
        axes[1, cloud].set_xlabel(r"$t/(10^7 M)$")
    axes[0, 0].set_yscale("log")
    axes[0, 1].ticklabel_format(axis="y", style="plain", useOffset=False)
    for axis in axes.flat:
        axis.set_xlim(0, times[-1]/1e7)
        axis.grid(alpha=0.2)
    figure.suptitle("Two-cloud outer-amplitude evolution", fontsize=17, y=0.98)
    seed_labels = []
    for seed in seeds:
        mantissa, exponent = ("%.2e" % seed).split("e")
        coefficient = "" if float(mantissa) == 1 else r"%g\times " % float(mantissa)
        seed_labels.append(coefficient + r"10^{%d}" % int(exponent))
    figure.text(0.5, 0.925,
                r"Fixed Kerr background: $\alpha=%.2f$, $\chi=%.2f$; " % (reference["alpha"], reference["spin"]) +
                r"$A_{211}(0)=%s$, $A_{322}(0)=%s$" % tuple(seed_labels), ha="center", fontsize=11)
    figure.legend([handles[k] for k in ("linear", "omiya", "cosine", "chiral")],
                  [labels[k] for k in ("linear", "omiya", "cosine", "chiral")],
                  loc="upper center", bbox_to_anchor=(0.5, 0.897), ncol=4, frameon=False)
    joint_domain = reference["spatial_settings"]
    omiya_domain = reports["omiya"]["spatial_settings"]
    note = ("Shared outer-amplitude seeds. Joint domain [%g, %g]; Omiya [%g, %g] in r*/M." %
            (joint_domain["inner"], joint_domain["outer"], omiya_domain["inner"], omiya_domain["outer"]))
    partial = ["%s: %d/%d steps" % (labels[k], m["steps"], len(times)-1)
               for k, m in models.items() if not m["complete"]]
    figure.text(0.5, 0.075, note, ha="center", fontsize=9, color="#444444")
    figure.text(0.5, 0.04, "Preview — " + "; ".join(partial) if partial else
                "Full common time interval; exploratory finite-domain comparison.",
                ha="center", fontsize=9, color="#444444")
    for extension in ("png", "pdf"):
        path = str(output) + "." + extension
        figure.savefig(path, dpi=180)
        print(path)
    plt.close(figure)
    np.savetxt(str(output) + "-trajectories.csv", np.column_stack(columns),
               delimiter=",", header=",".join(names), comments="")
    metadata = {"linear_model": "Independent exponential amplitude growth; no interactions or BH evolution",
                "equation": "dA_i/d(t/M)=Im(omega_i)*A_i; energy exponent is twice the amplitude exponent",
                "linear_rates": rates.tolist(), "normalization": "shared initial outer amplitudes, not matched occupations",
                "alpha": reference["alpha"], "spin": reference["spin"], "models": summary,
                "input_sha256": inputs, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(str(output) + "-comparison.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for model in ("cosine", "chiral", "omiya"):
        parser.add_argument("--" + model, required=True, type=Path, help="Trajectory output prefix")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    plot_growth_comparison(args.cosine, args.chiral, args.omiya, args.output, args.allow_partial)
