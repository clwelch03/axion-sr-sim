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
#SBATCH --output=runs/adiabatic/preliminary_20260922/logs/slurm-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export AXION_RUN_SECONDS=13800
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion_mpl_${SLURM_JOB_ID}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. runs/adiabatic/preliminary_20260922/branch/continue_cosine.jl
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python runs/adiabatic/preliminary_20260922/branch/plot_branch.py
