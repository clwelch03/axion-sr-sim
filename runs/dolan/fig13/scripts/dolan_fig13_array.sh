#!/bin/bash
#SBATCH -J dolan_td_arr
#SBATCH -A pc_heptheory
#SBATCH -p lr6
#SBATCH -q lr_normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -o runs/dolan/fig13/logs/%x_%A_%a.out
#SBATCH -e runs/dolan/fig13/logs/%x_%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: submit with --array, e.g. sbatch --array=0-14 $0" >&2
    exit 1
fi

# grid resolution: default 8, override at submit time with --export=ALL,GRID_N=4
grid_spacing_n="${GRID_N:-8}"

alpha=$(printf '%.2f' "$(echo "0.30 + 0.01 * ${SLURM_ARRAY_TASK_ID}" | bc -l)")
alpha_tag=${alpha/./} # 0.42 -> 042, for the filename
output="runs/dolan/fig13/sweep/dolan_fig13_a099_alpha${alpha_tag}_n${grid_spacing_n}_t1000000.npz"
mkdir -p "$(dirname "${output}")"
echo "array task ${SLURM_ARRAY_TASK_ID}: alpha=${alpha} -> ${output}"

module load miniconda3/22.11.1-gcc-11.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BHsuper          # numpy/scipy/mpmath/matplotlib

PYTHONPATH=src python -u -m dolan.measure_gamma \
    --output "${output}" \
    --time-final 1000000.0 \
    --spin 0.99 --alpha "${alpha}" --m 1 --l-max 5 --grid-spacing-n "${grid_spacing_n}" \
    --x-left -200.0 --x-right 511.0 \
    --x-pml-center -100.0 --pml-width 10.0 \
    --probe-x 20.0 --sample-every-steps 10 \
    --checkpoint-interval 10000.0

echo