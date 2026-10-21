#!/bin/bash
#SBATCH --job-name=axion_double_primary
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/doubling-%j.out
set -euo pipefail
snapshot_directory=$1
output_prefix=$2
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion-doubling-matplotlib-${SLURM_JOB_ID}"
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
/global/home/users/clwelch/.juliaup/bin/julia --project=. -e \
    'using Test; include("src/adiabatic/multiphase/load.jl"); include("tests/adiabatic/multiphase/test_evolution.jl"); include("tests/adiabatic/multiphase/test_amplitude_target.jl")'
/global/home/users/clwelch/.juliaup/bin/julia --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/double_primary.jl \
    inputs/initial.jls "$output_prefix"
PYTHONPATH=src /global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
    runs/adiabatic/multiphase_two_cloud/scripts/plot_doubling.py "$output_prefix" \
    --steps "${AXION_DOUBLING_STEPS:-64}"
