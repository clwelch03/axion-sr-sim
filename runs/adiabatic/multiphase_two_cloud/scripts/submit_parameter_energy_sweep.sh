#!/bin/bash
#SBATCH --job-name=axion_parameters
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/parameters-%A_%a.out
set -euo pipefail
[[ $# -eq 2 || $# -eq 3 ]] || { echo 'usage: SNAPSHOT OUTPUT_DIRECTORY [--validate]' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
output_directory=$(realpath -m "$2")
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$snapshot_directory"
sha256sum --quiet -c source-manifest.sha256
if [[ ${3:-} == --validate ]]; then
    /global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
        --startup-file=no --compiled-modules=existing --project=. \
        runs/adiabatic/multiphase_two_cloud/scripts/check_parameter_energy_sweep.jl \
        "$snapshot_directory/inputs" "$output_directory"
else
    [[ $# -eq 2 ]] || { echo 'unknown option' >&2; exit 2; }
    /global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
        --startup-file=no --compiled-modules=existing --project=. \
        runs/adiabatic/multiphase_two_cloud/scripts/parameter_energy_target.jl \
        "$snapshot_directory/inputs/initial.jls" "$snapshot_directory/inputs/sweep.toml" \
        "${SLURM_ARRAY_TASK_ID:?Submit with --array=1-N%4}" "$output_directory"
fi
