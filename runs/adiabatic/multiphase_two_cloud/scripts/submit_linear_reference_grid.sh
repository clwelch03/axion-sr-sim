#!/bin/bash
#SBATCH --job-name=axion_linear_grid
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/linear-grid-%A_%a.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'usage: SNAPSHOT OUTPUT_DIRECTORY' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
output_directory=$(realpath -m "$2")
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$snapshot_directory"
sha256sum --quiet -c source-manifest.sha256
/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
    --startup-file=no --compiled-modules=existing --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/linear_reference_grid.jl \
    "$snapshot_directory/inputs/linear_grid.toml" \
    "${SLURM_ARRAY_TASK_ID:?Submit with --array=1-12%4}" "$output_directory"
