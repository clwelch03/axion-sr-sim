"""Catherine's default matplotlib plot style for this project. Import it before making a figure:

    import matplotlib_style   # with scripts/ on the path, or: sys.path.insert(0, "scripts")

Importing sets the rcParams below (nothing else to call). Settings, from her own plotting configuration:
Computer Modern serif text, 2 pt lines, a six-colour cycle (cornflowerblue, forestgreen, maroon, goldenrod, firebrick,
mediumorchid), large axis labels (30 pt) and tick labels (24 pt), ticks on all four sides pointing inward with minor
ticks, a framed semi-transparent legend, 7 x 7 inch figures, tight saved bounding box.

Two differences from the configuration this was taken from:
  - LaTeX text rendering (`text.usetex`) is switched on only when a `latex` executable is on PATH. Without one,
    matplotlib would fail at the first label. On Lawrencium, load a TeX Live module first (the old hard-coded
    texlive/2016 PATH entry belongs to the retired sl-7 software tree and no longer exists).
  - Where the original set a value twice (tick label size 22, then 24) the later value, 24, is used.

The Julia figures (scripts/plot_*.jl, CairoMakie) do not use this file; they have their own theme in each script.
"""
import shutil

import matplotlib as mpl
from cycler import cycler

label_fontsize = 30   # axis-label size, also exported for use in titles or text added by hand
latex_available = shutil.which("latex") is not None

mpl.rcParams.update({
    # text: Computer Modern. With usetex the LaTeX Computer Modern is used; without it, the CMU fonts if installed.
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "font.sans-serif": ["CMU Sans Serif", "DejaVu Sans", "Bitstream Vera Sans", "Lucida Grande", "Verdana", "Geneva",
                        "Lucid", "Arial", "Helvetica", "Avant Garde", "sans-serif"],
    "text.usetex": latex_available,

    # lines and colours
    "lines.linewidth": 2.0,
    "lines.antialiased": True,
    "lines.dashed_pattern": (2.8, 1.5),
    "lines.dashdot_pattern": (4.8, 1.5, 0.8, 1.5),
    "lines.dotted_pattern": (1.1, 1.1),
    "lines.scale_dashes": True,
    "axes.prop_cycle": cycler("color", ["cornflowerblue", "forestgreen", "maroon", "goldenrod", "firebrick",
                                        "mediumorchid"]),

    # axes and ticks: ticks on all four sides, inward, with minor ticks
    "axes.linewidth": 1.0,
    "axes.labelsize": label_fontsize,
    "axes.labelpad": 10.0,
    "xtick.top": True, "ytick.right": True,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.size": 10, "ytick.major.size": 10,
    "xtick.minor.size": 5, "ytick.minor.size": 5,
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.minor.width": 0.75, "ytick.minor.width": 0.75,
    "xtick.major.pad": 8, "ytick.major.pad": 8,
    "xtick.labelsize": 24, "ytick.labelsize": 24,

    # legend
    "legend.fontsize": 18,
    "legend.frameon": True,
    "legend.framealpha": 0.8,
    "legend.fancybox": True,
    "legend.borderpad": 0.4,         # whitespace inside the legend border
    "legend.labelspacing": 0.5,      # vertical space between entries
    "legend.handlelength": 1.5,      # length of the legend lines
    "legend.handleheight": 0.7,
    "legend.handletextpad": 0.5,     # space between a legend line and its text
    "legend.borderaxespad": 0.5,     # space between the axes and the legend edge
    "legend.columnspacing": 2.0,

    # figure and saving
    "figure.figsize": (7, 7),
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})
