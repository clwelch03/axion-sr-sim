#!/bin/bash
# Submit with --dependency=afterany:<current-job-id> if another run is active.
#SBATCH --job-name=axion_branch_long
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=runs/adiabatic/fig2_20260922/logs/slurm-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export AXION_RUN_SECONDS=13800
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion_mpl_${SLURM_JOB_ID}"
export AXION_TARGET_AMPLITUDE="${AXION_TARGET_AMPLITUDE:-13}"
export AXION_PREDICTOR="${AXION_PREDICTOR:-cubic}" AXION_HISTORY_SPACING="${AXION_HISTORY_SPACING:-0.025}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. src/adiabatic/drivers/continue_branch.jl \
    cosine runs/adiabatic/fig2_20260922/seeds/preliminary_cosine_20260922.jls \
    runs/adiabatic/fig2_20260922/results/cosine_branch
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
