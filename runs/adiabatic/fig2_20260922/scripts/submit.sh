#!/bin/bash
#SBATCH --job-name=axion_branch
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=runs/adiabatic/fig2_20260922/logs/slurm-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export AXION_TARGET_AMPLITUDE="${AXION_TARGET_AMPLITUDE:-13}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. src/adiabatic/single_cloud/drivers/continue_branch.jl \
    cosine "${AXION_SEED_PATH:-runs/adiabatic/fig2_20260922/seeds/preliminary_cosine_20260922.jls}" \
    "${AXION_OUTPUT_PREFIX:-runs/adiabatic/fig2_20260922/results/cosine_branch}"
