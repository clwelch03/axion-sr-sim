#!/bin/bash
#SBATCH --job-name=axion_qcd_pilot
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/qcd-pilot-%A_%a.out
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
    runs/adiabatic/multiphase_two_cloud/scripts/qcd_profile_pilot.jl \
    "$snapshot_directory/inputs/pilots.toml" "${SLURM_ARRAY_TASK_ID:?Submit --array=1-2%2}" \
    "$output_directory"
