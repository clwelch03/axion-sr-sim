# Dolan Figure 13 time-domain runs

## Start here

This campaign measures linear quasi-bound growth rates from time-domain probe
signals for comparison with Dolan Fig. 13 and the eigenvalue/Leaver benchmarks.
For the standard comparison, use the `sweep/*_n6_t1000000.npz` records
($M/6$ spacing, the plotter's default). The `n4` and `n8` records are resolution
comparisons; `peak/dolan_fig13_a099_alpha042_t100000.npz` is the earlier,
shorter peak run. These NPZ files are local scratch artifacts, absent from a
fresh Git checkout.

From the repository root in the Python environment, plot with available
time-domain records or analyze one record:

```bash
PYTHONPATH=src python -m dolan.fig13_growth_rates --grid-spacing-n 6
PYTHONPATH=src python -m dolan.analyze_timedomain_run runs/dolan/fig13/sweep/dolan_fig13_a099_alpha042_n6_t1000000.npz
```

The plot command also computes eigenvalue/Leaver comparisons and writes
`plots/fig13_growth_rates.png`. The recipes below regenerate single points
or the sweep; select `GRID_N=6` for the standard comparison. To resume, run
the underlying `dolan.measure_gamma` command with the same numerical settings
and output path plus `--resume`. The submission scripts do not add that flag
automatically, so resubmitting one starts a fresh calculation.

**Validation limit:** a fixed final time does not give uniform first-overtone
accuracy across alpha, especially near low-growth regions and the threshold.
The [reproduction record](../../../docs/dolan-repro.md) contains dated
measurements; [HANDOFF.md](../../../HANDOFF.md) owns project status.

## Layout and preservation

`peak/` holds single-point results, `sweep/` holds the alpha/resolution sweep,
`scripts/` holds submission recipes, and `logs/` holds Slurm output.
Run commands from the repository root.

Heavy time-domain evolutions and their outputs. Per the repo `.gitignore`,
everything regenerable lives here but stays *out* of git: the run `.npz`, the
SLURM `.out`/`.err`, and any `.png`. Submit scripts (`.sh`) and this README
are tracked — they record how to regenerate the outputs.

## Dolan Fig. 13 time-domain point (a=0.99, alpha=0.42)

The growth-rate-peak point of Dolan 1212.1477 Fig. 13, evolved to t=10^5 M.
Massive-field Gaussian initial data (Eq. 45), lowest channel only; the fixed-x
(x=20M) complex probe series of every angular channel is saved for the
frequency-filter Gamma extraction.

- **Output:** `peak/dolan_fig13_a099_alpha042_t100000.npz` (gitignored, on scratch)
- **Cost:** ~1h07m single-core on Lawrencium `lr6` (job 25713497)
- **Analyse with:** `PYTHONPATH=src python -m dolan.analyze_timedomain_run <npz>` (raw spectrum + Gamma table)

### Submit script (`scripts/dolan_fig13_timedomain.sh`)

```bash
mkdir -p runs/dolan/fig13/logs runs/dolan/fig13/peak
sbatch runs/dolan/fig13/scripts/dolan_fig13_timedomain.sh
# Longer single-point run:
sbatch runs/dolan/fig13/scripts/dolan_fig13_timedomain_t1e6.sh
```

`measure_gamma.py` checkpoints every 5000 M, so a wall-time kill can be resumed
with the same command plus `--resume`.

## Full Fig. 13 time-domain curve (a=0.99, alpha sweep)

The whole curve as a SLURM array, one alpha per task -- script
`scripts/dolan_fig13_array.sh` (tracked alongside this recipe).
Each task is an independent t=10^6 M run into `runs/dolan/fig13/sweep/`.

- **alpha map:** `alpha = 0.30 + 0.01 * SLURM_ARRAY_TASK_ID` (task 0 -> 0.30,
  task 12 -> 0.42 peak, task 14 -> 0.44), matching `fig13_growth_rates.py`.
- **grid:** the `GRID_N` env var drives BOTH `--grid-spacing-n` and the output
  filename tag `_n${GRID_N}`. Submit an n=4 set (fast, ~8h/task) or n=8 set
  (accurate, ~31h/task):

  ```bash
  mkdir -p runs/dolan/fig13/logs runs/dolan/fig13/sweep
  sbatch --export=ALL,GRID_N=4 --array=0-14 runs/dolan/fig13/scripts/dolan_fig13_array.sh
  sbatch --export=ALL,GRID_N=8 --array=0-14 runs/dolan/fig13/scripts/dolan_fig13_array.sh
  ```

- **Why n=4 is enough:** the eigenvalue convergence test (`solve_qbs` vs grid n)
  gives only 0.6% discretization error at the alpha=0.42 peak for n=4, well below
  the ~3% time-domain extraction floor (filter ringing + n=1 SNR). Cost ~ n^2, so
  n=4 is 1/4 of n=8 for effectively the same curve. Pick grid n by that cheap
  eigenvalue test, not by defaulting high.
- **Analyse:** `PYTHONPATH=src python -m dolan.analyze_timedomain_run runs/dolan/fig13/sweep/<one>.npz` per point.
- **Caveat:** a fixed 10^6 M does not give uniform n=1 accuracy across alpha --
  the growth rate (hence SNR) shrinks at low alpha and near the n=1 zero-crossing
  (alpha~0.44), so the ends are noisier than the peak.
