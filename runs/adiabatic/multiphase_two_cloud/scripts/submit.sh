#!/bin/bash
#SBATCH --job-name=axion_joint_response
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/response-%j.out
set -euo pipefail
snapshot_directory=$1
output_prefix=$2
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
if [[ ${AXION_RUN_TESTS:-1} == 1 ]]; then
    /global/home/users/clwelch/.juliaup/bin/julia --project=. tests/adiabatic/runtests.jl
fi
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/response_and_evolution.jl "$output_prefix"
