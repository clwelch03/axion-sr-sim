#!/bin/bash
#SBATCH --job-name=axion_bvp_speed
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:50:00
#SBATCH --output=runs/adiabatic/radial_bvp/logs/speed-%j.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo "usage: submit_speed.sh VARIANT OUTPUT_PREFIX" >&2; exit 2; }
export RADIAL_BVP_ROOT="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp"
snapshot_directory="$RADIAL_BVP_ROOT/work/job-$SLURM_JOB_ID"
mkdir "$snapshot_directory"
cp "$RADIAL_BVP_ROOT"/scripts/*.jl "$snapshot_directory/"
cp "$RADIAL_BVP_ROOT/Project.toml" "$RADIAL_BVP_ROOT/Manifest.toml" "$snapshot_directory/"
bash "$RADIAL_BVP_ROOT/scripts/run_julia.sh" --compiled-modules=existing \
    "$snapshot_directory/speed_benchmark.jl" "$@"
