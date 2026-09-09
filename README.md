# axion-sr-sim

From-scratch reproduction of Dolan 1212.1477 Sec. III: the coupled-ℓ tower for a
massive scalar quasi-bound state on Kerr, evolved in the time domain, with an
eigenvalue cross-check and a Leaver continued-fraction benchmark.

- **`docs/dolan-repro.md`** — living status doc: where the reproduction stands,
  what has been validated, what is next.
- **`docs/DESIGN.md`** — design doc: why the code is shaped the way it is.
- **`CLAUDE.md`** — working conventions and how to run the solver, tests, and
  `runs/`.

## Quickstart

```
conda activate BHsuper                      # numpy, scipy, mpmath, matplotlib
python src/fig13_growth_rates.py            # eigenvalue Fig. 13 (fast)
python src/measure_gamma.py --help          # time-domain evolution driver
```

Modules live flat in `src/`; run scripts as `python src/<name>.py` from the repo
root. Heavy time-domain runs go through SLURM — see `runs/README.md`.
