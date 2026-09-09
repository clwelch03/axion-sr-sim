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
    --spin 0.99 --alpha 0.42 --m 1 --l-max 5 --n 6 \
    --x-left -200.0 --x-right 511.0 \
    --x-pml-center -100.0 --pml-width 10.0 \
    --probe-x 20.0 --sample-every-steps 10 \
    --checkpoint-interval 5000.0
```

`measure_gamma.py` checkpoints every 5000 M, so a wall-time kill can be resumed
with the same command plus `--resume`.
