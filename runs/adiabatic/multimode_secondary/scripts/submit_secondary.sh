#!/bin/bash
# Submit from the repository root. One writer per output prefix; rerunning resumes.
#SBATCH --job-name=axion_secondary
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multimode_secondary/logs/secondary-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    src/adiabatic/drivers/continue_secondary.jl \
    "${AXION_PRIMARY_PATH:-runs/adiabatic/fig2_20260922/results/cosine_branch.jls}" \
    "${AXION_OUTPUT_PREFIX:-runs/adiabatic/multimode_secondary/results/cosine_secondary}"
