#!/bin/bash
#SBATCH -J dolan_td_1e6
#SBATCH -A pc_heptheory
#SBATCH -p lr6
#SBATCH -q lr_normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#
# Longer Dolan Fig. 13 time-domain point (a=0.99, alpha=0.42), t=10^6 M.
# 10x the t=10^5 run (~1h07m) => ~11h wall. Purpose: confirm the n=1 overtone
# growth-rate deficit (ratio 0.87 at 1e5 M) is run-length/SNR-limited -- 10x
# longer accumulates ~10x the growth (0.011 -> ~0.11 in ln|A|^2 for n=1), which
# should dominate the ~0.01 filter-ringing wiggle and drive the ratio toward 1.0.
# Analyse with analyze_timedomain_run.py on the resulting .npz.
# Checkpoints every 5000 M (resumable with --resume + same command).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

module load miniconda3/22.11.1-gcc-11.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BHsuper          # numpy/scipy/mpmath/matplotlib

python -u src/measure_gamma.py \
    --output runs/dolan_fig13_a099_alpha042_t1000000.npz \
    --time-final 1000000.0 \
    --spin 0.99 --alpha 0.42 --m 1 --l-max 5 --n 6 \
    --x-left -200.0 --x-right 511.0 \
    --x-pml-center -100.0 --pml-width 10.0 \
    --probe-x 20.0 --sample-every-steps 10 \
    --checkpoint-interval 5000.0

echo "done."
