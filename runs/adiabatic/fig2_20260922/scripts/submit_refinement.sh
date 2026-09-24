#!/bin/bash
#SBATCH --job-name=axion_refine
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=runs/adiabatic/fig2_20260922/logs/refinement-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
unset AXION_SEED_PATH
if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then
    export AXION_POTENTIAL=quartic AXION_DIFFERENCE_STEP=1e-10
else
    export AXION_POTENTIAL=sextic AXION_DIFFERENCE_STEP=1e-8
fi
export AXION_NEWTON_MAX_ITERATIONS=12
export AXION_SAMPLE_SPACING="${AXION_SAMPLE_SPACING:-0.1}"
export AXION_RUN_SECONDS=13800
/global/home/users/clwelch/.juliaup/bin/julia --project=. runs/adiabatic/fig2_20260922/scripts/refine_polynomial.jl
