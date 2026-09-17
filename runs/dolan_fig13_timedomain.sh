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
#
# Dolan Fig. 13 time-domain point at the a=0.99, alpha=0.42 growth-rate peak.
# Evolves the massive-field Gaussian initial data (Eq. 45) to t=10^5 M and
# saves the fixed-x (x=20M) complex probe series of every angular channel.
# Analyse the resulting .npz with analyze_timedomain_run.py (raw spectrum +
# frequency-filter Gamma).  ~1h07m single-core on lr6 (job 25713497).
#
# Reproduces runs/dolan_fig13_a099_alpha042_t100000.npz (gitignored: *.npz).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# radial1d needs only numpy/scipy/mpmath/matplotlib -- the BHsuper conda env has them.
module load miniconda3/22.11.1-gcc-11.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BHsuper

PYTHONPATH=src python -u -m dolan.measure_gamma \
    --output runs/dolan_fig13_a099_alpha042_t100000.npz \
    --time-final 100000.0 \
    --spin 0.99 --alpha 0.42 --m 1 --l-max 5 --n 6 \
    --x-left -200.0 --x-right 511.0 \
    --x-pml-center -100.0 --pml-width 10.0 \
    --probe-x 20.0 --sample-every-steps 10 \
    --checkpoint-interval 5000.0

echo "done."
