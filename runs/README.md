# Simulation campaigns

Run submission commands from the repository root. Each campaign keeps its
recipes, results, seeds, and logs together. Reusable solvers remain in `src/`;
figures remain in `plots/`.

Each campaign's **Start here** section identifies its current inputs, results,
validation limits, and reproduction commands. For project-wide status and
priorities, use [HANDOFF.md](../HANDOFF.md).

| Campaign | Contents and reproduction |
|---|---|
| [Dolan Fig. 13](dolan/fig13/README.md) | Single-point time-domain runs and the alpha/resolution sweep |
| [Adiabatic Fig. 2, 2026-09-22](adiabatic/fig2_20260922/README.md) | Cosine, quartic, and sextic continuation, refinement, and diagnostics |
| [Small-amplitude checks](adiabatic/small_amplitude/README.md) | Initial continuation, profiles, and nonlinear convergence results |
| [Amplitude evolution](adiabatic/amplitude_evolution/README.md) | Merged 426-point derivative table, interpolation checks, and single-cloud A0(t) |
| [Multimode secondary](adiabatic/multimode_secondary/README.md) | 2211.01949 Sec. III: secondary cloud, driven fluxes, and Figs. 6–7 at Mμ = 0.42 and 0.3 |
| [Background continuation](adiabatic/background_continuation/README.md) | Local primary, secondary, and driven-field solves in alpha and spin, with route/return and radial-step checks |
| [QCD spin-down](adiabatic/qcd_spindown/README.md) | Preliminary coupled 20-solar-mass BH/211/322 evolution at m_a=1e-12 eV, compared with linear nonrelativistic growth |
| [QCD 10x10 grid](adiabatic/qcd_grid_20261012/README.md) | Alpha <= 0.45, 10-million-year cap; stopped before disconnect after a low-alpha pilot failure; frozen inputs and restart notes |

## Conventions

- `scripts/`: campaign drivers, submission scripts, and plotting recipes.
- `results/`: numerical tables and checkpoints; Dolan uses `peak/` and `sweep/`.
- `seeds/`: preserved inputs for restarts and comparisons.
- `diagnostics/`: experiments grouped by purpose, with their scripts and outputs.
- `logs/`: Slurm output, named with job IDs. Create this directory before submission.
- `README.md`: current inputs/results, validation limits, and reproduction or
  restart commands at the top; dated experiments below or in `HISTORY.md`.

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
