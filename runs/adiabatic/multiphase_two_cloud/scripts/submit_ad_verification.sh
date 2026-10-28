#!/bin/bash
#SBATCH --job-name=axion_ad_verify
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/ad-verify-%j.out
set -euo pipefail
[[ $# -eq 3 || $# -eq 4 ]] || { echo 'usage: SNAPSHOT STAGE OUTPUT_PREFIX [SEED_JLS]' >&2; exit 2; }
snapshot_directory=$1
stage=$2
output_prefix=$3
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. \
    src/adiabatic/multiphase/autodiff_vs_optimized_test.jl \
    "$stage" "$snapshot_directory/inputs/initial.jls" "$output_prefix" "${@:4}"
