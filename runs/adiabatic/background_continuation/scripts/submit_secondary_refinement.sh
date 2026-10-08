#!/bin/bash
#SBATCH --job-name=axion_secondary_refinement
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=runs/adiabatic/background_continuation/logs/secondary-refine-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/background_continuation/scripts/refine_secondary_backgrounds.jl
