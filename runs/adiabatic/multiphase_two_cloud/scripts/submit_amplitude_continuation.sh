#!/bin/bash
#SBATCH --job-name=axion_amplitude
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/amplitude-%j.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'usage: SNAPSHOT CAMPAIGN_DIRECTORY' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
campaign_directory=$(realpath "$2")
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1

# Verify and use the frozen code and seed; write results to the live campaign.
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/continue_amplitude.jl \
    "$campaign_directory" "$snapshot_directory/inputs/initial.jls"
