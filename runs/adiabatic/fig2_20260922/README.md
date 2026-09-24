# Adiabatic Figure 2 campaign — 2026-09-22

Preliminary cosine, quartic, and sextic branches at spin=0.99 and alpha=0.42.
These saved branches do not establish radial/angular discretization convergence.
See [HISTORY.md](HISTORY.md) for experiments, validation, and past job IDs;
those entries describe settings at the time, not necessarily today's defaults.

## Files

- `scripts/`: continuation, refinement, targeted quartic solves, and plotting.
- `results/`: accepted `*_branch.*`, additional `*_refined.*`, and quartic gap points.
- `seeds/`: initial cosine checkpoint and frozen diagnostic inputs.
- `diagnostics/`: extrapolation, speed, and high-amplitude comparisons.
- `logs/`: Slurm output. Figures go to `plots/preliminary_fig2_20260922.{png,pdf}`.

Snapshot at the layout migration (2026-09-24); CSVs remain authoritative:

| Result table | Saved points | Maximum amplitude |
|---|---:|---:|
| `cosine_branch.csv` | 178 | 11.9851915 |
| `quartic_branch.csv` | 126 | 6 |
| `quartic_gap.csv` | 11 | 5.4 |
| `quartic_refined.csv` | 235 | 6 |
| `sextic_branch.csv` | 87 | 10 |
| `sextic_refined.csv` | 154 | 9.10906286 |

## Reproduce or resume

From the repository root, create the logs directory before submitting:

```bash
mkdir -p runs/adiabatic/fig2_20260922/logs
sbatch runs/adiabatic/fig2_20260922/scripts/submit.sh
# Alternatively, a four-hour cosine continuation followed by plotting:
sbatch runs/adiabatic/fig2_20260922/scripts/submit_long.sh
# Independent polynomial branch (choose quartic or sextic):
sbatch --export=ALL,AXION_POTENTIAL=quartic runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
# Refine the two polynomial branches after their writers finish:
sbatch runs/adiabatic/fig2_20260922/scripts/submit_refinement.sh
```

The two cosine commands are alternatives. Never run two writers against the
same output prefix. Refinement should read completed production checkpoints.
Both continuation drivers resume their `results/*_branch.jls` when present;
otherwise they use `seeds/preliminary_cosine_20260922.jls`. Polynomial startup
re-solves the first seed point for the selected potential. `AXION_SEED_PATH`
can select an explicit seed. Cosine also accepts `AXION_OUTPUT_PREFIX`.

Current driver defaults: cosine uses ForwardDiff Jacobians, polynomial uses
finite differences; both use the secant predictor and a 12-correction limit.
`AXION_JACOBIAN`, `AXION_PREDICTOR`, `AXION_HISTORY_SPACING`,
`AXION_DIFFERENCE_STEP`, and `AXION_NEWTON_MAX_ITERATIONS` select alternatives.
`AXION_TARGET_AMPLITUDE` and `AXION_RUN_SECONDS` set continuation bounds;
submission scripts may set their own defaults. Existing serialized solutions
retain their saved physical/grid settings. Check the selected script and record
environment overrides before launching a new calculation.

Plot the existing data without running a solver:

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

The plot merges branch, refined, and gap CSVs, keeping production values at
shared amplitudes. Use the BHsuper Python environment. Julia recipes use the
repository `Project.toml` and `Manifest.toml` (`julia --project=.`).
