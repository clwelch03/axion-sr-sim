#!/bin/bash
#SBATCH --job-name=axion_solver_checks
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=runs/adiabatic/radialsolvertest/logs/checks-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=src/adiabatic/radialsolvertest src/adiabatic/radialsolvertest/check_methods.jl
