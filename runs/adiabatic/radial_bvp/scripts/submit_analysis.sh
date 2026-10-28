#!/bin/bash
#SBATCH --job-name=axion_bvp_check
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=runs/adiabatic/radial_bvp/logs/analysis-%j.out
set -euo pipefail
export RADIAL_BVP_ROOT="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp"
snapshot_directory="$RADIAL_BVP_ROOT/work/job-$SLURM_JOB_ID"
mkdir "$snapshot_directory"
cp "$RADIAL_BVP_ROOT"/scripts/*.jl "$snapshot_directory/"
cp "$RADIAL_BVP_ROOT/Project.toml" "$RADIAL_BVP_ROOT/Manifest.toml" "$snapshot_directory/"
bash "$RADIAL_BVP_ROOT/scripts/run_julia.sh" "$snapshot_directory/compare_results.jl" "$@"
