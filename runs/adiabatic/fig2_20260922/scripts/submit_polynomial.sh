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
export AXION_POTENTIAL="${AXION_POTENTIAL:-quartic}"
# Preserve this recipe's previous finite-difference and history-spacing defaults.
export AXION_JACOBIAN="${AXION_JACOBIAN:-finite}"
export AXION_HISTORY_SPACING="${AXION_HISTORY_SPACING:-0.05}"
/global/home/users/clwelch/.juliaup/bin/julia --project=. src/adiabatic/single_cloud/drivers/continue_branch.jl \
    "$AXION_POTENTIAL" "${AXION_SEED_PATH:-runs/adiabatic/fig2_20260922/seeds/preliminary_cosine_20260922.jls}" \
    "${AXION_OUTPUT_PREFIX:-runs/adiabatic/fig2_20260922/results/${AXION_POTENTIAL}_branch}"
