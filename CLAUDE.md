# CLAUDE.md — axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver benchmark. Extracted from the GRChombo fork
(`clwelch03/GRChombo-Axion-SR`, `Examples/AxionSuperradiance/radial1d/`); GRChombo
is no longer part of this line of work.

The repo now has two lines of work: the Dolan time-domain tower (Python,
`src/dolan/`) and the adiabatic nonlinear-cloud method of Omiya, Takahashi &
Tanaka (Julia, `src/adiabatic/`), which is the active one. The multimode
extension (Omiya et al. 2211.01949, Sec. III) is at the design stage.

Start with `docs/README.md` for the documentation index and `HANDOFF.md` for
where to resume (check its date; it can lag the docs). Current adiabatic state:
`docs/adiabatic-repro.md` (checkpoints) and `docs/adiabatic-branch-work-2026-09-28.md`
(latest session: dense derivatives, cosine flux zero, open items).
`docs/multimode-adiabatic-design.md` is a proposal only; no multimode solver
exists. `docs/dolan-repro.md` contains Dolan entry points and dated validation
history; `docs/DESIGN.md` is its unfinished design draft. Dated docs and reviews
(`*-review*.md`, `nonlinear-*.md`) are records, not current defaults. This file
is about how to work in the repo.

Reference PDFs are in `docs/papers/`. Cite equation numbers from those PDFs,
not from arXiv HTML (the numbers can differ).

## Working style on this project

- **Cheap decisive test over thorough survey.** Before proposing a mechanism for
  a discrepancy, ask what single measurement would falsify it. Prefer tests with
  an exact analytic answer (4th-order stencil ratio, Thomas vs `np.linalg.solve`,
  single-channel `rhs` vs a plain `b/D(r)` divide, spheroidal Leaver).
- **Falsify before mechanism.** Don't build the explanation before ruling out the
  boring cause. Example already in the history: the "~1% short of Leaver" was not
  a missing tower contribution — it was benchmarking against flat `Λ=l(l+1)`
  instead of spheroidal Leaver. One runtime comparison settled it.
- **Check against saved data, and say how far a check reaches.** State what a
  check does and does not establish (interpolation accuracy on a fixed grid is
  not radial/angular convergence; a local solve check is not discretization
  accuracy). Label preliminary results as such.
- **`evolve_metric=0 ⇒ time-independent flux`** is the canonical falsification
  principle carried over from the fork: if a knob that should not affect an
  observable does affect it, the observable is contaminated (this is how the
  PML-contaminated time-domain slope fit was caught — it drifted with PML width).

## Explain notation

If you introduce notation the paper or the repo does not already use (a symbol,
a combined phase, an abbreviation, a label for a quantity), define it the first
time it appears: say that it is your notation, give the definition, and say what
each piece means, including which of two similar symbols (for example the
Boyer–Lindquist azimuth versus the field φ₁) is meant. Do not assume the user
will recognize it. Example of what to avoid: writing "φ₁ depends on t and φ only
through ψ₁" with ψ₁ undefined. Write instead: "define ψ₁ ≡ ω₁t − m₁φ (my
notation, not the paper's), the combination in the exponent of the primary
ansatz (53), where φ is the Boyer–Lindquist azimuth." When a paper's notation
and the code's names differ (e.g. paper `m_n` versus code `m` / `fundamental_m`),
say which is which.

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
  resume a killed run with the same command plus `--resume`. Use `--help` for
  the full argument set and `runs/dolan/fig13/README.md` for submission recipes.
- **Analyse a finished time-domain run:**
  `PYTHONPATH=src python -m dolan.analyze_timedomain_run <run.npz>` — raw probe spectrum plus the
  frequency-filter Γ table against the runtime Leaver benchmark.
- **Leaver benchmark alone:** `PYTHONPATH=src python -m common.leaver --help` (continued fraction,
  mpmath; use `spheroidal=True` for the apples-to-apples comparison).

Heavy runs go through SLURM — see `runs/README.md` for the submit scripts.

### Julia (adiabatic)

Run from the repo root with `julia --project=. <script>`. Library code loads
through `src/adiabatic/load.jl`. Heavy runs go through Slurm, with Newton
residual tolerance and spatial settings recorded in checkpoints. Only one writer
per checkpoint prefix; an existing output is always resumed.

- **Continue a stationary branch (all potentials):**
  `AXION_TARGET_AMPLITUDE=13 julia --project=. src/adiabatic/drivers/continue_branch.jl cosine|quartic|sextic seed.jls output_prefix`.
  `AXION_TARGET_AMPLITUDE` is required. `continue_cosine.jl` and
  `continue_polynomial.jl` are superseded and pending removal; do not use them
  for new work.
- **Energy derivatives along a branch:** `src/adiabatic/drivers/branch_energy_derivatives.jl`
  (sharded by Slurm array; merge with `runs/adiabatic/amplitude_evolution/scripts/merge_derivatives.jl`).
- **Amplitude evolution:** `src/adiabatic/amplitude_evolution/` (`interpolation.jl`
  Hermite interpolation of the branch table, `time_integration.jl`,
  `branch_fitting.jl`, `branch_sensitivity.jl`); recipe
  `runs/adiabatic/amplitude_evolution/scripts/evolve_cosine.jl`.

Open caveats as of 2026-09-28: radial, outer-radius, mode-truncation, and angular
convergence at large amplitude are unchecked, so the cosine flux zero near
A0 ≈ 12.96 is not yet interpreted physically. Energy derivatives above
A0 = 11.985 needed a follow-up extension run. The seven `fig2_20260922/diagnostics/`
scripts run unchanged only at commit `aa17cb5`.

### Tests

`tests/dolan/test_hydrogenic.py` contains analytic normalization, overtone, threshold,
broadcasting, and input checks for the small-alpha benchmark. Run with
`PYTHONPATH=src python -m unittest discover -s tests/dolan` (also pytest-compatible).
Julia adiabatic checks: `julia --project=. tests/adiabatic/runtests.jl`
(`test_*.jl` are unit tests; `check_*.jl` are longer physical checks, some
needing Slurm).
The validation checks described in `solver.py`
(4th-order stencil convergence, Thomas vs dense solve, single-channel `rhs`
identity, tower ℓ-leakage) were run ad hoc; a grid-convergence script
(`convergence_test.py`) lives on the old fork branch
`feat/eigenmode-seed-and-backreaction`. Turning those into a real `tests/`
suite remains a pending task.

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

### Adiabatic modules

Load `src/adiabatic/load.jl` for the library. `angular.jl` and `tortoise.jl`
provide the shared basis and geometry. `stationary/` contains potentials,
radial equations, boundaries, integration, matching, and continuation predictors;
`observables/energy.jl` integrates energy. `amplitude_evolution/` interpolates the
saved branch table and integrates the slow amplitude A0(t). `drivers/` contains small- and
finite-amplitude continuation (`continue_branch.jl` for all three potentials),
branch energy derivatives, profile inspection, and the Julia Fig. 13 adapter.
Campaign submit scripts, results, and flux diagnostics live under `runs/adiabatic/`.
See `docs/adiabatic-design.md` for the complete map and conventions.

## runs/ layout

`runs/` groups simulations by campaign: `dolan/fig13/`,
`adiabatic/fig2_20260922/`, `adiabatic/small_amplitude/`, and
`adiabatic/amplitude_evolution/` (dense derivative tables, seeds, evolution).
See `runs/README.md` for the index, migration map, and preservation policy.

Campaign `scripts/` directories contain recipes; `results/` (or Dolan's
`peak/` and `sweep/`) contain numerical outputs; `seeds/` holds preserved inputs;
`logs/` holds Slurm output. Existing CSV/JLS reference results remain tracked.
Large NPZ files, logs, figures, temporary writes, and new `work/` scratch
experiments are ignored. Keep campaign READMEs and recipes current.
