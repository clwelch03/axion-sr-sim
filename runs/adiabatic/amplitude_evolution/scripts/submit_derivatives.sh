#!/bin/bash
# Submit from the repository root after freezing the input checkpoint.
# One writer per output shard; rerunning resumes accepted derivative points.
#SBATCH --job-name=axion_derivatives
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3%4
#SBATCH --time=01:00:00
#SBATCH --output=runs/adiabatic/amplitude_evolution/logs/derivatives-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    src/adiabatic/single_cloud/drivers/branch_energy_derivatives.jl \
    runs/adiabatic/amplitude_evolution/seeds/cosine_branch_20260924.jls \
    "runs/adiabatic/amplitude_evolution/results/cosine_derivatives_part${SLURM_ARRAY_TASK_ID}" \
    "${SLURM_ARRAY_TASK_ID}" 4
