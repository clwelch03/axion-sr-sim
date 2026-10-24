#!/bin/bash
#SBATCH --job-name=axion_direct_profile
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/direct-profile-%j.out
set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch submit_direct_profile.sh SNAPSHOT_DIRECTORY OUTPUT_PREFIX" >&2
    exit 2
fi
snapshot_directory=$1
output_prefix=$2
if [[ ! -d "$snapshot_directory" ]]; then
    printf 'Snapshot directory does not exist: <%s>\n' "$snapshot_directory" >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/direct_profile.jl \
    "$snapshot_directory/inputs" "$output_prefix"
