# Zoom on the spin bump of the a0 = 0.9, alpha0 = 0.3, f_a = 5.7e17 GeV cosine run and its continuation.
# Joins the chain of runs (original -> t1e11 restart -> later restart) and plots, against linear time,
# the spin with the |211> threshold a = 4 Omega / (1 + 4 Omega^2), Omega = Re omega_211 (current-mass units),
# and the two occupations. The |322> threshold, the same formula with Omega = Re omega_322 / 2 (azimuthal number 2),
# lies near a = 0.55, far below the bump, so the spin panel has a broken vertical axis: the zoom on the bump above,
# a narrow strip around the |322> threshold below.
#
# A finished segment is read from results/<prefix>_diagnostics.csv. A segment that is still running has no CSV
# (the driver writes it at the end), so it is read from its Slurm log instead: each "t = ..." line is one rhs
# evaluation with t, A, b, j. The log has no Re omega_211 and no occupations, so for a log segment
#   - each threshold is the straight-line fit of that threshold against mass b over the bump in the finished
#     segments, evaluated at the log's b (an extrapolation in b; drawn dotted to mark it),
#   - the occupations are N_i = N_i(last CSV row) * (A_i / A_i(last CSV row))^2 (weak-field scaling; drawn dotted).
#
# Run from the package root:
#   python scripts/plot_spin_bump.py <segment> [<segment> ...] [--output plots/name.png] [--start-time 1e10]
# where each segment is a results/ prefix (CSV) or a path to a Slurm log ending in .out.
import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("segments", nargs="+")
parser.add_argument("--output", default="plots/spin_bump_a0.9_alpha0.3_fa0.0466873.png")
parser.add_argument("--start-time", type=float, default=1e10)
arguments = parser.parse_args()


def accepted_rows(times):
    """Rows that are accepted states: the driver logs every rhs evaluation, and an accepted state is
    re-evaluated at the same time as the row before it (the filter of scripts/plot_spindown.jl)."""
    earliest_later_time = np.append(np.minimum.accumulate(times[::-1])[::-1][1:], np.inf)
    keep = (times <= earliest_later_time) & (np.append(np.nan, times[:-1]) == times)
    keep[0] = True
    return keep


def read_csv_segment(prefix):
    table = np.genfromtxt(f"results/{prefix}_diagnostics.csv", delimiter=",", names=True)
    table = table[accepted_rows(table["time"])]
    frequency_211, half_frequency_322 = table["Re_omega1"], table["Re_omega2"] / 2
    return dict(time=table["time"], spin=table["a"], mass=table["b"],
                threshold=4 * frequency_211 / (1 + 4 * frequency_211**2),
                threshold_322=4 * half_frequency_322 / (1 + 4 * half_frequency_322**2),
                occupation_211=table["N1"], occupation_322=table["N2"],
                amplitude_211=table["A1"], amplitude_322=table["A2"], from_log=False)


log_line = re.compile(r"^t = (\S+)\s+A = \[(\S+), (\S+)\].*?b = (\S+)\s+j = (\S+)")


def read_log_segment(path):
    rows = [[float(value) for value in match.groups()]
            for match in map(log_line.match, open(path)) if match]
    rows = np.array(rows)
    rows = rows[accepted_rows(rows[:, 0])]
    return dict(time=rows[:, 0], amplitude_211=rows[:, 1], amplitude_322=rows[:, 2], mass=rows[:, 3],
                spin=rows[:, 4] / rows[:, 3]**2, from_log=True)


segments = [read_log_segment(name) if name.endswith(".out") else read_csv_segment(name)
            for name in arguments.segments]
finished = [segment for segment in segments if not segment["from_log"]]
joined = {key: np.concatenate([segment[key] for segment in finished]) for key in finished[0] if key != "from_log"}

# threshold against mass over the bump, for log segments
in_bump = joined["time"] > 3e10
threshold_fits = {}
for name in ("threshold", "threshold_322"):
    threshold_fits[name] = np.polyfit(joined["mass"][in_bump], joined[name][in_bump], 1)
    fit_error = np.max(np.abs(np.polyval(threshold_fits[name], joined["mass"][in_bump]) - joined[name][in_bump]))
    print(f"{name}(b) straight-line fit over t > 3e10: max error {fit_error:.1e}")
for segment in segments:
    if segment["from_log"]:
        for name in threshold_fits:
            segment[name] = np.polyval(threshold_fits[name], segment["mass"])
        segment["occupation_211"] = joined["occupation_211"][-1] * (segment["amplitude_211"] / joined["amplitude_211"][-1])**2
        segment["occupation_322"] = joined["occupation_322"][-1] * (segment["amplitude_322"] / joined["amplitude_322"][-1])**2

figure = plt.figure(figsize=(8, 8.5))
grid = figure.add_gridspec(3, 1, height_ratios=[3, 0.7, 3], hspace=0.08)
spin_axis = figure.add_subplot(grid[0])
lower_spin_axis = figure.add_subplot(grid[1], sharex=spin_axis)   # strip around the |322> threshold
occupation_axis = figure.add_subplot(grid[2], sharex=spin_axis)
for segment in segments:
    shown = segment["time"] >= arguments.start_time
    time = segment["time"][shown] / 1e10
    derived_style = ":" if segment["from_log"] else "-"
    spin_axis.plot(time, segment["spin"][shown], color="black", lw=1.8)
    spin_axis.plot(time, segment["threshold"][shown], color="tab:green", lw=1.2, ls=derived_style)
    lower_spin_axis.plot(time, segment["threshold_322"][shown], color="tab:blue", lw=1.2, ls=derived_style)
    lower_spin_axis.plot(time, segment["spin"][shown], color="black", lw=1.8)   # enters the strip only near the threshold
    occupation_axis.plot(time, segment["occupation_211"][shown], color="tab:green", ls=derived_style)
    occupation_axis.plot(time, segment["occupation_322"][shown], color="tab:blue", ls=derived_style)
    print(f"segment to t = {segment['time'][-1]:.4e}: a = {segment['spin'][-1]:.6f}, threshold = {segment['threshold'][-1]:.6f}, "
          f"322 threshold = {segment['threshold_322'][-1]:.6f}, "
          f"N = ({segment['occupation_211'][-1]:.3g}, {segment['occupation_322'][-1]:.4g}), b = {segment['mass'][-1]:.6f}")
plateau_spin = joined["spin"][np.searchsorted(joined["time"], 1e10)]
spin_axis.axhline(plateau_spin, color="gray", lw=0.8, ls="--")
spin_axis.plot([], [], color="black", lw=1.8, label="spin $a = J/M^2$")
spin_axis.plot([], [], color="tab:green", lw=1.2, label=r"$|211\rangle$ threshold ($\omega_{211} = \Omega_H$)")
spin_axis.plot([], [], color="gray", lw=0.8, ls="--", label=f"plateau value {plateau_spin:.4f}")
if any(segment["from_log"] for segment in segments):
    spin_axis.plot([], [], color="tab:green", ls=":", label="from the running job's log (see script header)")
spin_axis.plot([], [], color="tab:blue", lw=1.2, label=r"$|322\rangle$ threshold ($\omega_{322} = 2\Omega_H$), lower strip")
spin_axis.set_ylabel("black-hole spin")
spin_axis.legend(loc="best", fontsize=9)
all_threshold_322 = np.concatenate([segment["threshold_322"][segment["time"] >= arguments.start_time] for segment in segments])
lower_spin_axis.set_ylim(all_threshold_322.min() - 0.004, all_threshold_322.max() + 0.004)
# broken-axis marks between the zoom and the strip
spin_axis.spines["bottom"].set_visible(False)
lower_spin_axis.spines["top"].set_visible(False)
spin_axis.tick_params(bottom=False, labelbottom=False)
lower_spin_axis.tick_params(labelbottom=False)
break_marks = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=9, linestyle="none", color="black", clip_on=False)
spin_axis.plot([0, 1], [0, 0], transform=spin_axis.transAxes, **break_marks)
lower_spin_axis.plot([0, 1], [1, 1], transform=lower_spin_axis.transAxes, **break_marks)
occupation_axis.plot([], [], color="tab:green", label=r"$|211\rangle$")
occupation_axis.plot([], [], color="tab:blue", label=r"$|322\rangle$")
occupation_axis.set_yscale("log")
occupation_axis.set_ylim(1e-6, 1e2)
occupation_axis.set_ylabel(r"occupation $\mathcal{N}_i$  [$f_a^2 M_0^2$]")
occupation_axis.set_xlabel(r"$t / (10^{10} M_0)$")
occupation_axis.legend(loc="best", fontsize=9)
figure.suptitle(r"self-interacting (cosine), $a_0 = 0.9$, $\alpha_0 = 0.3$, $f_a = 5.7\times10^{17}$ GeV: spin bump and continuation", fontsize=11)
figure.subplots_adjust(top=0.94, bottom=0.07, left=0.12, right=0.97)
figure.savefig(arguments.output, dpi=160)
print("wrote", arguments.output)
