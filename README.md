# axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver continued-fraction benchmark. A second
reproduction, Omiya, Takahashi & Tanaka [2201.04382](https://arxiv.org/abs/2201.04382),
is being built incrementally in `src/adiabatic/`.

- **`docs/dolan-repro.md`** — living status doc: where the reproduction stands,
  what has been validated, what is next.
- **`docs/adiabatic-repro.md`** — adiabatic method and first validation checkpoint.
- **`docs/DESIGN.md`** — design doc: why the code is shaped the way it is.
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

Cheap analytic checks: `PYTHONPATH=src python -m unittest discover -s tests`.
