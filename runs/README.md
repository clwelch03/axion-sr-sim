# Simulation campaigns

Run submission commands from the repository root. Each campaign keeps its
recipes, results, seeds, and logs together. Reusable solvers remain in `src/`;
figures remain in `plots/`.

| Campaign | Contents and reproduction |
|---|---|
| [Dolan Fig. 13](dolan/fig13/README.md) | Single-point time-domain runs and the alpha/resolution sweep |
| [Adiabatic Fig. 2, 2026-09-22](adiabatic/fig2_20260922/README.md) | Cosine, quartic, and sextic continuation, refinement, and diagnostics |
| [Small-amplitude checks](adiabatic/small_amplitude/README.md) | Initial continuation, profiles, and nonlinear convergence results |

## Conventions

- `scripts/`: campaign drivers, submission scripts, and plotting recipes.
- `results/`: numerical tables and checkpoints; Dolan uses `peak/` and `sweep/`.
- `seeds/`: preserved inputs for restarts and comparisons.
- `diagnostics/`: experiments grouped by purpose, with their scripts and outputs.
- `logs/`: Slurm output, named with job IDs. Create this directory before submission.
- `README.md`: current entry points and settings; `HISTORY.md`: dated experiments.

Record the code revision, submission command (including environment overrides),
job IDs, input seed, and output paths when starting or resuming a campaign.
Only one writer may use a given checkpoint prefix at a time. A restart continues
an existing campaign; use a new directory for a distinct configuration study.

## Version control and preservation

Track recipes, configuration, summaries, and selected reference tables and seeds.
Existing tracked CSV/JLS results are retained in this reorganization. Large Dolan
NPZ files, Slurm logs, and generated figures remain ignored. Temporary writes and
new scratch experiments under a campaign's `work/` directory are ignored. Do not
blanket-ignore or remove existing checkpoints: some are expensive restart inputs.

## Layout migration, 2026-09-24

The Slurm queue was empty before migration, including the dependent Figure 2
refresh. All 194 moved data/log artifacts were verified byte-for-byte with SHA-256.
No numerical results were deleted or regenerated during the move.

| Previous location | Current location |
|---|---|
| `runs/dolan_*.sh` | `runs/dolan/fig13/scripts/` |
| `runs/dolan_*.npz` | `runs/dolan/fig13/peak/` |
| `runs/fig13/*.npz` | `runs/dolan/fig13/sweep/` |
| Dolan `.out` / `.err` files | `runs/dolan/fig13/logs/` |
| `runs/adiabatic/preliminary_20260922/branch/` | `fig2_20260922/scripts/` and `fig2_20260922/results/` under `runs/adiabatic/` |
| Other `preliminary_20260922/` contents | Corresponding folders under `runs/adiabatic/fig2_20260922/` |
| `runs/adiabatic/preliminary_cosine_20260922.*` | `runs/adiabatic/fig2_20260922/seeds/` |
| Loose adiabatic continuation/profile/convergence files | `runs/adiabatic/small_amplitude/results/` |

Serialized provenance strings and historical logs retain their original paths.
Use this table to locate those inputs; current scripts use the new locations.
The old campaign README is preserved as `fig2_20260922/HISTORY.md`.
