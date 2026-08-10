#!/bin/bash
#SBATCH --job-name=towerG
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr4
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --array=0-9
#SBATCH --output=/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/nonlinear/gamma_sweep/towerG_%A_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/nonlinear/gamma_sweep/towerG_%A_%a.err
#
# Gamma(A, l_max) multi-mode growth sweep at the superradiant growth peak
# (a=0.99, alpha=0.42). One array task per (l_max, amplitude, lambda) point.
# Pure-numpy single-process RK4 -- no BLAS parallelism to exploit.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SWEEP=/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/nonlinear/gamma_sweep
GWENV=/global/scratch/projects/pc_heptheory/clwelch/software/envs/gwenv/bin/python

# --- sweep grid (index by SLURM_ARRAY_TASK_ID) ---------------------------------
#   task 0     : lambda=0 LINEAR baseline (Gamma_lin freqfilter anchor for differencing)
#   task 1-3   : l_max=1 CONTROL  -- Hermitian single-mode cubic => Gamma must stay ~+1.5e-7
#   task 4-6   : l_max=3          -- first cascade channel; Gamma(A) should start to fall
#   task 7-9   : l_max=5          -- l_max convergence of the growth suppression
LMAX=(1    1    1    1    3    3    3    5    5    5   )
AMP=( 2    1    2    4    1    2    4    1    2    4   )
LAM=( 0    0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02)

i=$SLURM_ARRAY_TASK_ID
L=${LMAX[$i]}; A=${AMP[$i]}; LM=${LAM[$i]}
OUT=$SWEEP/tower_gamma_lmax${L}_A${A}_lam${LM}.npz

echo "task $i: l_max=$L amp=$A lambda=$LM -> $OUT"
cd "$SWEEP"
$GWENV "$SWEEP/tower_gamma_run.py" --l_max "$L" --amp "$A" --lam "$LM" --out "$OUT"
