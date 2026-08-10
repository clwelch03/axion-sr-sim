#!/usr/bin/env python
"""Animate the 1+1D radial psi(x,t): the 1D analogue of the project-level
../../viz/cloud_movie.py's equatorial-slice movie (note: DIFFERENT directory
from this file's own radial1d/viz/ -- that one is the 3+1D GRChombo movie
tool, this one is 1D-radial-specific) -- Re[psi(x)] and |psi(x)| as a
function of the tortoise coordinate, one frame per recorded snapshot -> GIF.

Shows directly what the static-heatmap version couldn't as clearly: the
ingoing wavefront reaching the horizon, the bound-state envelope settling in
at the cloud peak, and (if PML reflection is present) any standing-wave
ripple in |psi| near x_match.

Run (needs a --save with --nsnap>0 from evolve.py; evolve.py lives one
directory up, in radial1d/ itself, so run this from radial1d/, not from
radial1d/viz/ -- or invoke evolve.py with its full path):
    cd ..   # into radial1d/
    $GWENV evolve.py --alpha 0.50 --a 0.90 --T 3000 --N 4000 --nsnap 300 \
        --save viz/pml_a050_field.npz
    cd viz
    $GWENV spacetime.py pml_a050_field.npz --out pml_a050_movie
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# palette (dataviz reference, light mode) -- matches ../viz/cloud_movie.py
SURFACE = "#fcfcfb"; INK = "#0b0b0b"; INK2 = "#52514e"; MUTED = "#898781"
RE_COLOR = "#2f6fb0"; AMP_COLOR = "#b0402f"; PML_FILL = "#e7e5e0"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz", help="output of evolve.py --save ... --nsnap N")
    ap.add_argument("--out", default=None, help="output basename (default: npz stem)")
    ap.add_argument("--fps", type=float, default=15.0)
    ap.add_argument("--stride", type=int, default=1, help="use every Nth snapshot")
    args = ap.parse_args()

    d = np.load(args.npz)
    if "field_psi" not in d:
        raise SystemExit(f"{args.npz} has no field snapshots -- rerun evolve.py with --nsnap > 0")

    t = d["field_t"][::args.stride]; x = d["x"]; psi = d["field_psi"][::args.stride]
    rp = float(d["rp"]); x_match = float(d["x_match"])
    out = args.out or args.npz.rsplit(".", 1)[0] + "_movie"

    amp = np.abs(psi)
    ymax = 1.05 * np.max(amp)   # fixed scale across all frames (decay/growth stays visible)

    frames = []
    for i in range(len(t)):
        fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=120, facecolor=SURFACE)
        ax.axvspan(x[0], x_match, color=PML_FILL, zorder=0, label="PML region")
        ax.axvline(x_match, color=MUTED, lw=0.9, ls="--")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.plot(x, psi[i].real, color=RE_COLOR, lw=1.2, label=r"Re $\psi$")
        ax.plot(x, amp[i], color=AMP_COLOR, lw=1.1, label=r"$|\psi|$")
        ax.plot(x, -amp[i], color=AMP_COLOR, lw=1.1)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(-ymax, ymax)
        ax.set_xlabel("tortoise x / M", fontsize=9, color=INK2)
        ax.set_ylabel(r"$\psi$", fontsize=9, color=INK2)
        ax.tick_params(colors=MUTED, labelsize=8)
        for s in ax.spines.values(): s.set_color(MUTED)
        ax.set_title("1+1D QBS evolution on Kerr (Dolan 2012 reduction)",
                     fontsize=10, color=INK, pad=8)
        ax.text(0.02, 0.96, f"t = {t[i]:7.1f} M", transform=ax.transAxes,
                fontsize=10, color=INK, va="top", ha="left",
                bbox=dict(fc=SURFACE, ec="none", alpha=0.85))
        ax.text(0.98, 0.04, rf"horizon ($r_+={rp:.3f}$) $\to$ off left edge",
                transform=ax.transAxes, fontsize=7.5, color=INK2, va="bottom",
                ha="right", bbox=dict(fc=SURFACE, ec="none", alpha=0.75))
        ax.legend(loc="upper right", fontsize=8, frameon=False, labelcolor=INK2)
        fig.tight_layout()
        fig.canvas.draw()
        img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        frames.append(Image.fromarray(img))
        plt.close(fig)
        if i % 20 == 0:
            print(f"  frame {i}/{len(t)} (t={t[i]:.0f})")

    gif = f"{out}.gif"
    frames[0].save(gif, save_all=True, append_images=frames[1:],
                   duration=int(1000 / args.fps), loop=0)
    print(f"wrote {gif} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
