#!/bin/bash
#SBATCH --job-name=axion_weak_secondary
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/weak-secondary-%j.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo "usage: submit_weak_secondary.sh SNAPSHOT OUTPUT_PREFIX" >&2; exit 2; }
snapshot_directory=$1
output_prefix=$2
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/weak_secondary_profile.jl \
    "$snapshot_directory/inputs/initial.jls" "$output_prefix"
