# Adiabatic Figure 2 campaign — 2026-09-22

## Start here

This campaign computes stationary cosine, quartic, and sextic branches at
spin 0.99 and alpha 0.42. Use `results/cosine_branch.{csv,jls}` for the
primary input to the multimode campaign. For the Figure 2 plot, use the
production `*_branch.csv` tables together with `quartic_gap.csv` and the
quartic/sextic `*_refined.csv` tables; the plotting script merges them and
retains production values at shared amplitudes.

From the repository root, redraw the saved tables without solving:

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

Output: `plots/preliminary_fig2_20260922.{png,pdf}`. To continue a branch,
use [Reproduce or resume](#reproduce-or-resume); existing output prefixes
resume automatically. `seeds/` supplies frozen inputs for fresh outputs and
diagnostic comparisons, rather than newer production results.

**Validation limit:** the branches are preliminary. The cosine primary is
known to fail radial-step convergence above approximately $A_0 = 9$, so the
flux zero near 12.96 is not an established physical saturation point. Full
high-amplitude boundary, angular, and mode-truncation convergence is also
unestablished. See [the handoff](../../../HANDOFF.md#known-caveats) for
project caveats and the [amplitude campaign](../amplitude_evolution/README.md)
for derivatives and single-cloud time evolution.

Older diagnostic scripts require the revision in
[Historical driver freeze](#historical-driver-freeze).
[HISTORY.md](HISTORY.md) and the
[September 28 work record](../../../docs/adiabatic-branch-work-2026-09-28.md)
retain their dated experiments and next-step lists.

## Files and current results

- `scripts/`: submission recipes, polynomial gap refinement, and plotting.
- `results/`: accepted `*_branch.*`, `*_refined.*`, and quartic gap points.
- `seeds/`: initial cosine checkpoint and frozen diagnostic inputs.
- `diagnostics/`: historical experiments; see the freeze point below.
- `logs/`: Slurm output. Figures go to `plots/preliminary_fig2_20260922.{png,pdf}`.

Saved tables on September 28 (CSVs remain authoritative):

| Result table | Saved points | Maximum amplitude |
|---|---:|---:|
| `cosine_branch.csv` | 232 | 13 |
| `quartic_branch.csv` | 126 | 6 |
| `quartic_gap.csv` | 11 | 5.4 |
| `quartic_refined.csv` | 235 | 6 |
| `sextic_branch.csv` | 87 | 10 |
| `sextic_refined.csv` | 154 | 9.10906286 |

Cosine job **26434700** added 54 points beyond A0=11.9851915. Energy remains
monotonic, and the signed net loss changes from negative to positive between
A0=12.959108 and 12.962707. The linear zero estimate is A0=12.96268 on this
spatial grid. Its physical interpretation requires spatial convergence checks.
Energy and flux derivatives live in the separate
[amplitude evolution campaign](../amplitude_evolution/README.md).

## Reproduce or resume

The unified driver is `src/adiabatic/drivers/continue_branch.jl`:

```bash
AXION_TARGET_AMPLITUDE=13 julia --project=. src/adiabatic/drivers/continue_branch.jl \
    cosine runs/adiabatic/fig2_20260922/seeds/preliminary_cosine_20260922.jls \
    runs/adiabatic/fig2_20260922/results/cosine_branch
```

The three positional arguments are potential (`cosine`, `quartic`, or `sextic`),
seed checkpoint, and output prefix. The target environment variable is required.
An existing output always resumes; the seed is used only for a fresh output.
A seed for another potential supplies its first point, which is re-solved.
Use a new output prefix to start over. Physical/grid settings come from the
selected checkpoint. Never run two writers against the same output prefix.

From the repository root:

```bash
mkdir -p runs/adiabatic/fig2_20260922/logs
sbatch runs/adiabatic/fig2_20260922/scripts/submit.sh
# Alternative: four-hour cosine continuation followed by plotting.
sbatch runs/adiabatic/fig2_20260922/scripts/submit_long.sh
# Independent polynomial branch; target should be chosen explicitly.
sbatch --export=ALL,AXION_POTENTIAL=quartic,AXION_TARGET_AMPLITUDE=6 runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
# Fill polynomial gaps after production writers finish.
sbatch runs/adiabatic/fig2_20260922/scripts/submit_refinement.sh
```

The unified driver defaults to ForwardDiff, a secant predictor, 12 Newton
corrections, tolerance 1e-9, and a 1000-second runtime. Every continuation
attempt obeys `delta A0 <= min(0.25, 0.05 A0)`, including resumed attempts.
`AXION_MAX_AMPLITUDE_STEP` and `AXION_MAX_FRACTIONAL_STEP` override these limits.
`AXION_JACOBIAN`, `AXION_PREDICTOR`, `AXION_HISTORY_SPACING`,
`AXION_DIFFERENCE_STEP`, `AXION_NEWTON_MAX_ITERATIONS`, `AXION_INITIAL_STEP`,
and `AXION_RUN_SECONDS` set the corresponding solver/continuation controls.

Recipes retain their own defaults: short cosine targets 13; long cosine uses
cubic prediction with history spacing 0.025 and 13800 seconds; polynomial
uses finite differences, history spacing 0.05, and target 0.1 unless overridden.
The short cosine and polynomial recipes accept `AXION_SEED_PATH` and
`AXION_OUTPUT_PREFIX`. The long recipe uses the production paths above.
Polynomial refinement uses `load_branch_setup` and `branch_diagnostics` from
the unified driver but retains its independent midpoint gap-filling algorithm.
It writes `*_refined.*`, with `AXION_SAMPLE_SPACING=0.1` by default.

Plot without running a solver (BHsuper Python environment):

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

The plot merges branch, refined, and gap CSVs, retaining production values at
shared amplitudes. Julia recipes use `julia --project=.`.

## Historical driver freeze

The old cosine and polynomial drivers were replaced by `continue_branch.jl`.
Commit **aa17cb5** is the last revision at which all seven scripts under
`diagnostics/`, plus `scripts/solve_quartic_target.jl` and
`scripts/fill_quartic_gap.jl`, run unchanged. Reproduce those finished experiments
from a separate checkout of that revision; their old-driver references are
historical and are not current entry points.
