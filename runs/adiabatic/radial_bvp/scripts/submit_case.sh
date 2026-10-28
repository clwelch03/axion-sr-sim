#!/bin/bash
#SBATCH --job-name=axion_radial_bvp
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:55:00
#SBATCH --output=runs/adiabatic/radial_bvp/logs/bvp-%j.out
set -euo pipefail
[[ $# -ge 4 ]] || { echo "usage: submit_case.sh CASE TOL DT OUTPUT_PREFIX [RESTART]" >&2; exit 2; }
export RADIAL_BVP_ROOT="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp"
snapshot_directory="$RADIAL_BVP_ROOT/work/job-$SLURM_JOB_ID"
mkdir "$snapshot_directory"
cp "$RADIAL_BVP_ROOT"/scripts/*.jl "$snapshot_directory/"
cp "$RADIAL_BVP_ROOT/Project.toml" "$RADIAL_BVP_ROOT/Manifest.toml" "$snapshot_directory/"
bash "$RADIAL_BVP_ROOT/scripts/run_julia.sh" --compiled-modules=existing \
    "$snapshot_directory/solve_case.jl" "$@"
