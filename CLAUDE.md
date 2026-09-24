# CLAUDE.md — axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver benchmark. Extracted from the GRChombo fork
(`clwelch03/GRChombo-Axion-SR`, `Examples/AxionSuperradiance/radial1d/`); GRChombo
is no longer part of this line of work.

`docs/dolan-repro.md` is the living status doc for the reproduction — read it
first for where the physics stands. `docs/DESIGN.md` is the design doc (the
why). This file is about how to work in the repo.

## Working style on this project

- **Cheap decisive test over thorough survey.** Before proposing a mechanism for
  a discrepancy, ask what single measurement would falsify it. Prefer tests with
  an exact analytic answer (4th-order stencil ratio, Thomas vs `np.linalg.solve`,
  single-channel `rhs` vs a plain `b/D(r)` divide, spheroidal Leaver).
- **Falsify before mechanism.** Don't build the explanation before ruling out the
  boring cause. Example already in the history: the "~1% short of Leaver" was not
  a missing tower contribution — it was benchmarking against flat `Λ=l(l+1)`
  instead of spheroidal Leaver. One runtime comparison settled it.
- **`evolve_metric=0 ⇒ time-independent flux`** is the canonical falsification
  principle carried over from the fork: if a knob that should not affect an
  observable does affect it, the observable is contaminated (this is how the
  PML-contaminated time-domain slope fit was caught — it drifted with PML width).

## Naming style

Long informative variable names: `cloud_norm`, not `anorm`; `x_pml_center`, not
`xc`. Legacy short names in older code can stay; new code uses the long form.

## Say what you mean

Plain statements. Avoid mannered prose, overly flowerly language, and metaphor
for the sake of metaphor. Instead of "this point earns its keep", write "this point
still matters"; instead of "a dial worth turning", write "a parameter worth varying".

## Building from-scratch reproductions

This rule governs the whole repo. The user is a physicist and a relative novice
to numerical methods, and is using this reproduction to learn the method.

- Build incrementally, one piece at a time. The user follows and confirms each
  step.
- Do not hand back a finished module. Do not silently lift code. Do not skip the
  hard part to get a fast result.
- Each step gets a cheap decisive checkpoint before moving on (see above).

## How to run

Environment: conda env `BHsuper` (numpy, scipy, mpmath, matplotlib). On
Lawrencium: `module load miniconda3/22.11.1-gcc-11.4.0` then
`conda activate BHsuper`.

Packages live under `src/dolan/`, `src/adiabatic/`, and `src/common/`.
Use qualified imports (`from dolan.solver import Tower`) and run drivers from
the repo root with `PYTHONPATH=src python -m dolan.<module>`.
An editable install (`python -m pip install -e . --no-deps`) allows omitting
`PYTHONPATH=src`. Pytest includes `src` through `pyproject.toml`.
For adiabatic work, read `docs/adiabatic-design.md` for the numerical design and
`docs/adiabatic-repro.md` for current results and remaining checkpoints.

- **Eigenvalue growth rate (linear, fast):**
  `PYTHONPATH=src python -m dolan.fig13_growth_rates` — sweeps α at a=0.99, overlays `solve_qbs`
  eigenvalue points on spheroidal Leaver for n=0 and n=1. Writes a PNG.
- **Time-domain growth rate (the path that extends to nonlinear):**
  `PYTHONPATH=src python -m dolan.measure_gamma --output <run.npz> --time-final <T> --spin <a>
  --alpha <α> --m 1 --l-max 5 ...` — evolves Dolan's Gaussian initial data,
  saves the fixed-x complex probe series per channel. Checkpoints on an interval;
  resume a killed run with the same command plus `--resume`. Full argument set
  and the reference invocation are in `runs/README.md`.
- **Analyse a finished time-domain run:**
  `PYTHONPATH=src python -m dolan.analyze_timedomain_run <run.npz>` — raw probe spectrum plus the
  frequency-filter Γ table against the runtime Leaver benchmark.
- **Leaver benchmark alone:** `PYTHONPATH=src python -m common.leaver --help` (continued fraction,
  mpmath; use `spheroidal=True` for the apples-to-apples comparison).

Heavy runs go through SLURM — see `runs/README.md` for the submit scripts.

### Tests

`tests/dolan/test_hydrogenic.py` contains analytic normalization, overtone, threshold,
broadcasting, and input checks for the small-alpha benchmark. Run with
`PYTHONPATH=src python -m unittest discover -s tests/dolan` (also pytest-compatible).
Julia adiabatic checks: `julia --project=. tests/adiabatic/runtests.jl`.
The validation checks described in `solver.py`
(4th-order stencil convergence, Thomas vs dense solve, single-channel `rhs`
identity, tower ℓ-leakage) were run ad hoc; a grid-convergence script
(`convergence_test.py`) lives on the old fork branch
`feat/eigenmode-seed-and-backreaction`. Turning those into a real `tests/`
suite is a pending task (see `HANDOFF.md` / the repo-structure plan).

## Module map

Dolan modules are in `src/dolan/`; `leaver.py`, `hydrogenic.py`, and
`plotstyle.py` are in `src/common/`. Library (imported by drivers):
- `solver.py` — grid, Kerr tortoise `r(x)`, the `c^m_{jl}` coupling matrix,
  Eq. 25/26 coefficients, 4th-order stencils, Thomas solve, RK4/MoL. The `Tower`
  object is the core data model.
- `leaver.py` — Dolan-2007 continued-fraction QBS solver (benchmark).
- `qbs_eigenmode.py` — `solve_qbs`: the QBS as an eigenvector of the linear
  generator; gives the eigenvalue Γ.
- `frequency_filter.py` — Dolan's frequency filter (Eq. 46/49), pure importable.
- `plotstyle.py` — the repo default plot palette; sets the matplotlib colour
  cycle (five brand colours as C0..C4, topped up to ten distinct) on import.
- `fig13_timedomain_points.py` — sweeps the finished `measure_gamma.py` runs and
  extracts the frequency-filter Γ per α, selecting runs by the `grid_spacing_n`
  recorded inside each `.npz`. Imported by `fig13_growth_rates` for the overlay.

Drivers:
- `measure_gamma.py` — time-domain evolution + probe recording; imports
  `frequency_filter`.
- `fig13_growth_rates.py` — the Fig. 13 reproduction: Leaver curves, our
  eigenvalue points, and (via `fig13_timedomain_points`) the time-domain Γ
  points; `--grid-spacing-n 4 8` overlays several resolutions.
- `analyze_timedomain_run.py` — spectrum + Γ table from a saved probe record.

## runs/ layout

`runs/` groups simulations by campaign: `dolan/fig13/`,
`adiabatic/fig2_20260922/`, and `adiabatic/small_amplitude/`.
See `runs/README.md` for the index, migration map, and preservation policy.

Campaign `scripts/` directories contain recipes; `results/` (or Dolan's
`peak/` and `sweep/`) contain numerical outputs; `seeds/` holds preserved inputs;
`logs/` holds Slurm output. Existing CSV/JLS reference results remain tracked.
Large NPZ files, logs, figures, temporary writes, and new `work/` scratch
experiments are ignored. Keep campaign READMEs and recipes current.
