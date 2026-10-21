#!/bin/bash
#SBATCH --job-name=axion_omiya_compare
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/omiya-comparison-%j.out
set -euo pipefail
snapshot_directory=$1
joint_prefix=$2
output_prefix=$3
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/omiya_comparison.jl \
    "$joint_prefix" "$output_prefix"
