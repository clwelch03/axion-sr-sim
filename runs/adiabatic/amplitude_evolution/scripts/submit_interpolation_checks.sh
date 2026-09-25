#!/bin/bash
# Prepare queries with tests/adiabatic/check_branch_interpolation.jl first.
#SBATCH --job-name=axion_interpolation
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3%4
#SBATCH --time=00:30:00
#SBATCH --output=runs/adiabatic/amplitude_evolution/logs/interpolation-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    src/adiabatic/drivers/branch_energy_derivatives.jl \
    runs/adiabatic/amplitude_evolution/seeds/interpolation_checks_20260924.jls \
    "runs/adiabatic/amplitude_evolution/results/interpolation_checks_part${SLURM_ARRAY_TASK_ID}" \
    "${SLURM_ARRAY_TASK_ID}" 4
