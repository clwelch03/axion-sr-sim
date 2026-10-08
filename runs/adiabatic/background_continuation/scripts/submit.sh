#!/bin/bash
# Submit from the repository root. One job writes this campaign's checkpoints.
# Supply OUTPUT_DIRECTORY BRANCH_INDEX as sbatch script arguments.
#SBATCH --job-name=axion_background_check
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --output=runs/adiabatic/background_continuation/logs/nearby-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
[[ $# == 2 ]] || { echo "usage: sbatch submit.sh OUTPUT_DIRECTORY BRANCH_INDEX" >&2; exit 2; }
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/background_continuation/scripts/check_nearby_backgrounds.jl "$@"
