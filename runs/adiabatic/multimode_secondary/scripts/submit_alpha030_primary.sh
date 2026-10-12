#!/bin/bash
# Mμ = 0.3 cosine primary for the Fig. 7 / Fig. 6 (right) comparison.
# Run prepare_alpha030_seed.jl once first. One writer per output; rerunning resumes.
#SBATCH --job-name=axion_primary_a030
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=runs/adiabatic/multimode_secondary/logs/primary_a030-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
export AXION_TARGET_AMPLITUDE="${AXION_TARGET_AMPLITUDE:-0.06}"
export AXION_RUN_SECONDS="${AXION_RUN_SECONDS:-13800}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. src/adiabatic/single_cloud/drivers/continue_branch.jl \
    cosine runs/adiabatic/multimode_secondary/seeds/alpha030_linear_seed.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch
