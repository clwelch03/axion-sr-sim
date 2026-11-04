# axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver continued-fraction benchmark. A second
reproduction, Omiya, Takahashi & Tanaka [2201.04382](https://arxiv.org/abs/2201.04382),
uses stationary nonlinear profiles, continuation, energy, and flux
calculations for cosine, quartic, and sextic potentials in `src/adiabatic/`.
The multimode extension follows Omiya et al. 2211.01949 Sec. III, coupling a
secondary cloud and driven fields to the nonlinear primary on fixed Kerr.

Start with the [project handoff](HANDOFF.md) for current capabilities,
validation limits, and next steps. Use the [documentation index](docs/README.md)
for equations and research records, and the [campaign index](runs/README.md)
for saved results and reproduction commands.

- **[Dolan reproduction](docs/dolan-repro.md)** — current entry points and
  dated linear validation results.
- **[`docs/adiabatic-design.md`](docs/adiabatic-design.md)** — adiabatic equations,
  code structure, numerical choices, nonlinear matching, and observables.
- **[`docs/adiabatic-repro.md`](docs/adiabatic-repro.md)** — adiabatic validation
  records and observable conventions.
- **[Dolan design](docs/DESIGN.md)** — unfinished time-domain design draft.
- **[CLAUDE.md](CLAUDE.md)** — working conventions and how to run the solver, tests, and
  `runs/`.

## Setup

For the temporary Lawrencium-to-Perlmutter move, follow
[the transfer and setup guide](docs/perlmutter-handoff.md).

Run commands from the repository root. Python requires version 3.9 or newer
and NumPy, SciPy, mpmath, and Matplotlib. In an activated Python environment:

```bash
python -m pip install -e .
```

On Lawrencium the existing environment is `BHsuper`:

```bash
module load miniconda3/22.11.1-gcc-11.4.0
conda activate BHsuper
```

Adiabatic solves also require Julia and the repository project dependencies
(listed in `Project.toml`):

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Use `--project=.` for Julia commands. The checked-in `Manifest.toml` records
the dependency environment. Python-only Dolan calculations do not require Julia.

## Quickstart

```bash
PYTHONPATH=src python -m dolan.fig13_growth_rates            # eigenvalue Fig. 13 (fast)
PYTHONPATH=src python -m dolan.measure_gamma --help          # time-domain evolution driver
```

Packages live under `src/`:

- `dolan/`: time-domain solver, nonlinear extensions, eigenmode solver and drivers.
- `adiabatic/`: [model-based layout](src/adiabatic/README.md): shared `core/`,
  `single_cloud/`, `primary_secondary/`, and the planned `multiphase/` extension.
- `common/`: independent Leaver and hydrogenic benchmarks, plotting utilities.

Run the commands above from the repository root. Alternatively, install with
`python -m pip install -e . --no-deps` and omit `PYTHONPATH=src`.
Heavy time-domain runs go through SLURM — see `runs/README.md`.

The optional hydrogenic benchmark overlays the leading small-alpha rate:

```
PYTHONPATH=src python -m dolan.fig13_growth_rates --skip-eigenvalue --include-hydrogenic
```

This writes `plots/fig13_with_hydrogenic.png`. The approximation is extrapolated
over the Fig. 13 range; agreement near the peak is not assumed. The standalone
`common.hydrogenic.growth_rate(alpha, spin, radial_overtone=0, l=1, m=1)` returns signed
`M*omega_I`, the field-amplitude growth rate. It adapts the rate in
`clwelch03/axion-superradiance` (`src/superradiance.py`) with an explicit division
by two: [Ng et al., Eq. (1)](https://arxiv.org/html/1908.02312#S2.E1) gives the
occupation-number rate. Its radial overtone corresponds to Leaver's principal
number minus `l+1`. The horizon factor retains the source approximation
`omega_R = mu`, so its finite-alpha threshold is approximate too.

The adiabatic linear solver can be added to the same Fig. 13 axes as open diamonds:

```
PYTHONPATH=src python -m dolan.fig13_growth_rates --include-adiabatic
```

This requires the Julia project dependencies (`julia --project=.`); pass
`--julia /path/to/julia` if Julia is not on `PATH`. It solves both radial overtones
with angular modes `l=[1,3,5]`, using the Leaver real frequencies to select the
branches and starting the imaginary frequency at zero. The figure remains
`plots/fig13_growth_rates.png`; frequencies, residuals, and numerical settings
are saved in `plots/fig13_adiabatic_points.csv`. Add `--skip-eigenvalue` to omit
the ARPACK sweep, or combine with `--include-hydrogenic` for that comparison.
Time-domain points default to `M/6`; use `--grid-spacing-n` to select another resolution.

Redraw both preliminary Fig. 2 panels from the saved tables (no solve):

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

This writes `plots/preliminary_fig2_20260922.png` and `.pdf`. For continuation
and restart recipes, see the campaign README linked above.

Python checks: `PYTHONPATH=src python -m unittest discover -s tests/dolan`.

Julia adiabatic checks: `julia --project=. tests/adiabatic/runtests.jl`.
