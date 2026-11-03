#!/bin/bash
#SBATCH --job-name=axion_growth_energy
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/growth-energy-%j.out
set -euo pipefail
[[ $# -ge 3 && $# -le 4 ]] || { echo 'usage: COSINE_PREFIX CHIRAL_PREFIX OUTPUT_PREFIX [STRIDE]' >&2; exit 2; }
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
    --startup-file=no --compiled-modules=existing --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/export_growth_energies.jl "$@"
