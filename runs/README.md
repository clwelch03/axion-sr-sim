# radial1d time-domain runs

Heavy time-domain evolutions and their outputs. Per the repo `.gitignore`,
everything regenerable lives here but stays *out* of git: the run `.npz`, the
SLURM `.out`/`.err`, the `.sh` submit script, and any `.png`. Only this README
is committed — it records how to regenerate the rest.

## Dolan Fig. 13 time-domain point (a=0.99, alpha=0.42)

The growth-rate-peak point of Dolan 1212.1477 Fig. 13, evolved to t=10^5 M.
Massive-field Gaussian initial data (Eq. 45), lowest channel only; the fixed-x
(x=20M) complex probe series of every angular channel is saved for the
frequency-filter Gamma extraction.

- **Output:** `dolan_fig13_a099_alpha042_t100000.npz` (gitignored, on scratch)
- **Cost:** ~1h07m single-core on Lawrencium `lr6` (job 25713497)
- **Analyse with:** `analyze_timedomain_run.py <npz>` (raw spectrum + Gamma table)

### Submit script (`dolan_fig13_timedomain.sh`)

```bash
#!/bin/bash
#SBATCH -J dolan_td
#SBATCH -A pc_heptheory
#SBATCH -p lr6
#SBATCH -q lr_normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

module load miniconda3/22.11.1-gcc-11.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BHsuper          # numpy/scipy/mpmath/matplotlib

python -u measure_gamma.py \
    --output runs/dolan_fig13_a099_alpha042_t100000.npz \
    --time-final 100000.0 \
    --spin 0.99 --alpha 0.42 --m 1 --l-max 5 --grid-spacing-n 6 \
    --x-left -200.0 --x-right 511.0 \
    --x-pml-center -100.0 --pml-width 10.0 \
    --probe-x 20.0 --sample-every-steps 10 \
    --checkpoint-interval 5000.0
```

`measure_gamma.py` checkpoints every 5000 M, so a wall-time kill can be resumed
with the same command plus `--resume`.

## Full Fig. 13 time-domain curve (a=0.99, alpha sweep)

The whole curve as a SLURM array, one alpha per task -- script
`dolan_fig13_array.sh` (a `.sh`, so gitignored; kept here as the committed recipe).
Each task is an independent t=10^6 M run into `runs/fig13/`.

- **alpha map:** `alpha = 0.30 + 0.01 * SLURM_ARRAY_TASK_ID` (task 0 -> 0.30,
  task 12 -> 0.42 peak, task 14 -> 0.44), matching `fig13_growth_rates.py`.
- **grid:** the `GRID_N` env var drives BOTH `--grid-spacing-n` and the output
  filename tag `_n${GRID_N}`. Submit an n=4 set (fast, ~8h/task) or n=8 set
  (accurate, ~31h/task):

  ```bash
  mkdir -p runs/fig13                                   # for the SBATCH -o/-e logs
  sbatch --export=ALL,GRID_N=4 --array=0-14 runs/dolan_fig13_array.sh
  sbatch --export=ALL,GRID_N=8 --array=0-14 runs/dolan_fig13_array.sh
  ```

- **Why n=4 is enough:** the eigenvalue convergence test (`solve_qbs` vs grid n)
  gives only 0.6% discretization error at the alpha=0.42 peak for n=4, well below
  the ~3% time-domain extraction floor (filter ringing + n=1 SNR). Cost ~ n^2, so
  n=4 is 1/4 of n=8 for effectively the same curve. Pick grid n by that cheap
  eigenvalue test, not by defaulting high.
- **Analyse:** `analyze_timedomain_run.py runs/fig13/<one>.npz` per point.
- **Caveat:** a fixed 10^6 M does not give uniform n=1 accuracy across alpha --
  the growth rate (hence SNR) shrinks at low alpha and near the n=1 zero-crossing
  (alpha~0.44), so the ends are noisier than the peak.
