#!/bin/bash
#SBATCH --job-name=axion_references
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/reference-profiles-%j.out
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'usage: SNAPSHOT CAMPAIGN_DIRECTORY OUTPUT_PREFIX' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
campaign_directory=$(realpath "$2")
output_prefix=$(realpath -m "$3")
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1

# Run the reviewed source copy while reading preserved campaign inputs.
cd "$snapshot_directory"
sha256sum --quiet -c source-manifest.sha256
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/reference_profiles.jl \
    "$campaign_directory" "$output_prefix"
