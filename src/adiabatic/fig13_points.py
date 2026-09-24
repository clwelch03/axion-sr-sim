"""Run the Julia linear adiabatic solver for the two Fig. 13 radial branches."""

from pathlib import Path
import subprocess
import tempfile

import numpy as np


def growth_rate_points(alpha_values, references, spin, output_path, julia="julia"):
    """Return {label: complex frequencies}, also saving solver diagnostics to CSV.

    references is the branch-tracked Leaver map from the figure driver. Julia
    uses its real frequencies as initial guesses, but solves for growth itself.
    The retained angular modes are l=[1,3,5], m=1 for both radial overtones.
    """
    source_directory = Path(__file__).resolve().parent
    script = source_directory / "drivers" / "fig13_points.jl"
    project_root = source_directory.parents[1]
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels = ("n=0", "n=1")
    rows = []
    for overtone, label in enumerate(labels):
        for alpha in alpha_values:
            frequency = references[label][round(float(alpha), 3)]
            rows.append((overtone, alpha, frequency.real, frequency.imag))
    with tempfile.TemporaryDirectory(prefix="adiabatic-fig13-") as directory:
        seed_path = Path(directory) / "seeds.csv"
        np.savetxt(seed_path, rows, delimiter=",")
        subprocess.run(
            [julia, f"--project={project_root}", str(script), str(seed_path),
             str(output_path), str(spin), "0.05", "15"], check=True)

    points = np.genfromtxt(output_path, delimiter=",", names=True, ndmin=1)
    results = {}
    for overtone, label in enumerate(labels):
        branch = points[points["radial_overtone"] == overtone]
        if len(branch) != len(alpha_values) or not np.allclose(branch["alpha"], alpha_values):
            raise ValueError(f"Adiabatic output has missing or reordered points for {label}")
        frequencies = branch["omega_real"] + 1j * branch["omega_imag"]
        if not np.all(np.isfinite(frequencies)) or not np.all(branch["residual_norm"] <= 1e-10):
            raise ValueError(f"Adiabatic output contains unconverged points for {label}")
        results[label] = frequencies
    return results
