#!/bin/bash
#SBATCH --job-name=axion_profile_remake
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/profile-remake-%j.out
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'usage: SNAPSHOT OUTPUT_DIRECTORY check|joint|references' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
output_directory=$(realpath -m "$2")
run_mode=$3
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$snapshot_directory"
sha256sum --quiet -c source-manifest.sha256
julia_executable=/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia
if [[ $run_mode == check ]]; then
    "$julia_executable" --startup-file=no --compiled-modules=existing --project=. \
        tests/adiabatic/primary_secondary/test_comparison_secondary_energy.jl
    "$julia_executable" --startup-file=no --compiled-modules=existing --project=. \
        runs/adiabatic/multiphase_two_cloud/scripts/remake_profile_comparison.jl \
        "$snapshot_directory/inputs" "$output_directory" --check
elif [[ $run_mode == joint || $run_mode == references ]]; then
    driver=remake_profile_comparison.jl
    if [[ $run_mode == references ]]; then driver=remake_profile_references.jl; fi
    "$julia_executable" --startup-file=no --compiled-modules=existing --project=. \
        "runs/adiabatic/multiphase_two_cloud/scripts/$driver" \
        "$snapshot_directory/inputs" "$output_directory"
else
    echo 'unknown run mode' >&2
    exit 2
fi
