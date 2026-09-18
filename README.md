# axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver continued-fraction benchmark. A second
reproduction, Omiya, Takahashi & Tanaka [2201.04382](https://arxiv.org/abs/2201.04382),
is being built incrementally in `src/adiabatic/`.

- **`docs/dolan-repro.md`** — living status doc: where the reproduction stands,
  what has been validated, what is next.
- **[`docs/adiabatic-design.md`](docs/adiabatic-design.md)** — adiabatic equations,
  code structure, numerical choices, and nonlinear extension plan.
- **[`docs/adiabatic-repro.md`](docs/adiabatic-repro.md)** — adiabatic validation
  results, run commands, and remaining checkpoints.
- **`docs/DESIGN.md`** — Dolan time-domain design.
- **`CLAUDE.md`** — working conventions and how to run the solver, tests, and
  `runs/`.

## Quickstart

```
conda activate BHsuper                      # numpy, scipy, mpmath, matplotlib
PYTHONPATH=src python -m dolan.fig13_growth_rates            # eigenvalue Fig. 13 (fast)
PYTHONPATH=src python -m dolan.measure_gamma --help          # time-domain evolution driver
```

Packages live under `src/`:

- `dolan/`: time-domain solver, nonlinear extensions, eigenmode solver and drivers.
- `adiabatic/`: new reproduction of nonlinear profiles and slow cloud evolution.
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

Python checks: `PYTHONPATH=src python -m unittest discover -s tests/dolan`.

Julia adiabatic checks: `julia --project=. tests/adiabatic/runtests.jl`.
