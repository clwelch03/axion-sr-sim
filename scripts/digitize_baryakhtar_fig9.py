"""Digitize the |211> spin-down stage of Baryakhtar, Galanis, Lasenby & Simon, arXiv:2011.11646, Fig. 9 right panel
(10 M_sun, mu = 2.5e-12 eV) from docs/papers/2011.11646.pdf (version dated Feb 18, 2021, page 28).

Method: render page 28 at 600 dpi with Ghostscript, crop the right panel, calibrate the axes from the tick marks found
in the crop (y: rows of the 0.9 and 0.3 ticks; x: columns of the 1e2 and 1e10 ticks, log axis), and for each f take the
pixels whose colour is within COLOUR_TOLERANCE of that curve's legend colour, in a column window around its drop and
between a = 0.885 and a = 0.67 (so the shared a = 0.9 line, the threshold line, and the legend are excluded). Each
column gives one point: the mean row of the matching pixels. Reading error is about one line width, ~0.003 in a and
~0.02 in log10 t.

Only f = 1e11, 5e12, 5e14 and 5e15 GeV are digitized (the curves our plateau spin-down can be compared with).

Run from the package root (`module load python` on Perlmutter; needs numpy and Pillow, and gs on the path):
    python -I scripts/digitize_baryakhtar_fig9.py
Writes results/baryakhtar_fig9_right_digitized.csv with columns f_gev, log10_t_yr, spin.
"""
import os
import subprocess
import tempfile

import numpy as np
from PIL import Image

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
pdf_path = os.path.join(root, "docs", "papers", "2011.11646.pdf")
output_path = os.path.join(root, "results", "baryakhtar_fig9_right_digitized.csv")

with tempfile.TemporaryDirectory() as scratch:
    page_path = os.path.join(scratch, "page28.png")
    subprocess.run(["gs", "-q", "-sDEVICE=png16m", "-r600", "-dFirstPage=28", "-dLastPage=28", "-o", page_path, pdf_path],
                   check=True)
    page = np.asarray(Image.open(page_path).convert("RGB")).astype(int)
panel = page[400:1520, 2700:4600]          # right panel at 600 dpi (frame at columns 136-1847, rows 34-1092 of the crop)
dark = panel.sum(axis=2) < 250

# axis calibration from the tick marks
left_tick_rows = [row for row in range(40, 1080) if dark[row, 138:150].sum() >= 10 and dark[row, 160:200].sum() < 5]
bottom_tick_columns = [column for column in range(140, 1846)
                       if dark[1078:1090, column].sum() >= 8 and dark[1040:1070, column].sum() < 3]
assert len(left_tick_rows) == 4 and len(bottom_tick_columns) == 5, (left_tick_rows, bottom_tick_columns)
row_09, row_03 = left_tick_rows[0], left_tick_rows[-1]          # a = 0.9, 0.3
column_1e2, column_1e10 = bottom_tick_columns[0], bottom_tick_columns[-1]


def spin_of_row(row):
    return 0.9 - (row - row_09) * 0.6 / (row_03 - row_09)


def row_of_spin(spin):
    return row_09 + (0.9 - spin) * (row_03 - row_09) / 0.6


def log_time_of_column(column):
    return 2 + (column - column_1e2) * 8 / (column_1e10 - column_1e2)


# legend colours (median of the legend line segments) and the column window of each curve's drop
curves = {
    1e11: ((113, 113, 197), (1500, 1846)),
    5e12: ((100, 148, 236), (1050, 1500)),
    5e14: ((92, 171, 237), (520, 950)),
    5e15: ((187, 209, 237), (310, 650)),
}
COLOUR_TOLERANCE = 35
top_row, bottom_row = int(row_of_spin(0.885)), int(row_of_spin(0.67))

rows_out = []
for decay_constant, (colour, (column_start, column_stop)) in curves.items():
    distance = np.sqrt(((panel - np.array(colour)) ** 2).sum(axis=2))
    for column in range(column_start, column_stop):
        matches = np.where(distance[top_row:bottom_row, column] < COLOUR_TOLERANCE)[0]
        if len(matches):
            rows_out.append((decay_constant, log_time_of_column(column), spin_of_row(top_row + matches.mean())))

with open(output_path, "w") as output:
    output.write("f_gev,log10_t_yr,spin\n")
    for decay_constant, log_time, spin in rows_out:
        output.write(f"{decay_constant:.3e},{log_time:.5f},{spin:.5f}\n")
for decay_constant in curves:
    points = [(t, a) for f, t, a in rows_out if f == decay_constant]
    print(f"f = {decay_constant:.0e}: {len(points)} points, log10 t from {points[0][0]:.2f} to {points[-1][0]:.2f}, "
          f"a from {points[0][1]:.3f} to {points[-1][1]:.3f}")
print(f"calibration: a = 0.9 at row {row_09}, 0.3 at row {row_03}; t = 1e2 yr at column {column_1e2}, 1e10 at {column_1e10}")
print(f"wrote {output_path}")
