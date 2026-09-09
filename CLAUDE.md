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

Plain statements. Avoid flowery language for its own sake and metaphor for the
sake of metaphor. A design note should read like documentation, not an essay.

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

The modules live flat in `src/` and import each other by bare name
(`from solver import Tower`). Run them as `python src/<script>.py` from the repo
root — Python puts the script's own directory (`src/`) on the path, so the bare
imports resolve. (`pytest` finds them via `pythonpath = ["src"]` in
`pyproject.toml`; `pip install -e .` also works but is not needed.)

- **Eigenvalue growth rate (linear, fast):**
  `python src/fig13_growth_rates.py` — sweeps α at a=0.99, overlays `solve_qbs`
  eigenvalue points on spheroidal Leaver for n=0 and n=1. Writes a PNG.
- **Time-domain growth rate (the path that extends to nonlinear):**
  `python src/measure_gamma.py --output <run.npz> --time-final <T> --spin <a>
  --alpha <α> --m 1 --l-max 5 ...` — evolves Dolan's Gaussian initial data,
  saves the fixed-x complex probe series per channel. Checkpoints on an interval;
  resume a killed run with the same command plus `--resume`. Full argument set
  and the reference invocation are in `runs/README.md`.
- **Analyse a finished time-domain run:**
  `python src/analyze_timedomain_run.py <run.npz>` — raw probe spectrum plus the
  frequency-filter Γ table against the runtime Leaver benchmark.
- **Leaver benchmark alone:** `python src/leaver.py --help` (continued fraction,
  mpmath; use `spheroidal=True` for the apples-to-apples comparison).

Heavy runs go through SLURM — see `runs/README.md` for the submit scripts.

### Tests

There is no `tests/` suite yet. The validation checks described in `solver.py`
(4th-order stencil convergence, Thomas vs dense solve, single-channel `rhs`
identity, tower ℓ-leakage) were run ad hoc; a grid-convergence script
(`convergence_test.py`) lives on the old fork branch
`feat/eigenmode-seed-and-backreaction`. Turning those into a real `tests/`
directory is a pending task (see `HANDOFF.md` / the repo-structure plan).

## Module map

All modules are in `src/`. Library (imported by drivers):
- `solver.py` — grid, Kerr tortoise `r(x)`, the `c^m_{jl}` coupling matrix,
  Eq. 25/26 coefficients, 4th-order stencils, Thomas solve, RK4/MoL. The `Tower`
  object is the core data model.
- `leaver.py` — Dolan-2007 continued-fraction QBS solver (benchmark).
- `qbs_eigenmode.py` — `solve_qbs`: the QBS as an eigenvector of the linear
  generator; gives the eigenvalue Γ.
- `frequency_filter.py` — Dolan's frequency filter (Eq. 46/49), pure importable.

Drivers:
- `measure_gamma.py` — time-domain evolution + probe recording; imports
  `frequency_filter`.
- `fig13_growth_rates.py` — the eigenvalue Fig. 13 reproduction.
- `analyze_timedomain_run.py` — spectrum + Γ table from a saved probe record.

## runs/ layout

`runs/` holds heavy time-domain evolutions and their outputs.

- **Tracked:** `runs/README.md` (regeneration recipes) and the SLURM submit
  scripts (`*.sh`).
- **Gitignored (regenerable):** `*.npz` run outputs, `*.out`/`*.err` SLURM logs,
  `*.png`. `runs/fig13/` holds the per-α array outputs.

The `.sh` scripts plus `runs/README.md` are the definition of a reproduction —
keep them current when a run recipe changes.
