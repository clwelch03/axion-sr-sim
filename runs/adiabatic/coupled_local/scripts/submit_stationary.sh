#!/bin/bash
#SBATCH --job-name=axion_coupled_local
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --output=runs/adiabatic/coupled_local/logs/stationary-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/coupled_local/scripts/sample_stationary.jl "$@"
