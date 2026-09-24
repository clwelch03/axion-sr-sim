#!/bin/bash
#SBATCH --job-name=axion_solver_test
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-8%3
#SBATCH --output=runs/adiabatic/radialsolvertest/logs/case-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=src/adiabatic/radialsolvertest src/adiabatic/radialsolvertest/run_case.jl
