#!/bin/bash
#SBATCH --job-name=axion_polynomial
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=runs/adiabatic/fig2_20260922/logs/polynomial-%x-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export AXION_TARGET_AMPLITUDE="${AXION_TARGET_AMPLITUDE:-0.1}"
export AXION_RUN_SECONDS="${AXION_RUN_SECONDS:-1000}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. runs/adiabatic/fig2_20260922/scripts/continue_polynomial.jl
