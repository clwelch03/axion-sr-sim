#!/bin/bash
#SBATCH --job-name=qcd_grid
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-82%4
#SBATCH --output=runs/adiabatic/qcd_grid_20261012/logs/grid-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion-grid-matplotlib-${SLURM_JOB_ID}"
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
    runs/adiabatic/qcd_grid_20261012/scripts/campaign.py run "${SLURM_ARRAY_TASK_ID}"
