#!/bin/bash
# Submit from the repository root; one writer for the secondary output directory.
#SBATCH --job-name=axion_secondary_backgrounds
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --output=runs/adiabatic/background_continuation/logs/secondary-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/background_continuation/scripts/check_secondary_backgrounds.jl \
    runs/adiabatic/background_continuation/results/amplitude2052/nearby_backgrounds.jls \
    runs/adiabatic/background_continuation/results/amplitude2052/secondary
