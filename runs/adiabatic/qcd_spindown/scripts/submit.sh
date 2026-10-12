#!/bin/bash
#SBATCH --job-name=qcd_spindown
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/qcd_spindown/logs/stationary-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    src/adiabatic/drivers/qcd_stationary_point.jl "$@"
