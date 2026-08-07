#!/bin/bash
#SBATCH --job-name=stageE
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr4
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --output=stage_e_sweep_%j.out
#SBATCH --error=stage_e_sweep_%j.err
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
GWENV=/global/scratch/projects/pc_heptheory/clwelch/software/envs/gwenv/bin/python
cd /global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d
$GWENV stage_e_sweep.py
