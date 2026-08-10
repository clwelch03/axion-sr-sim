#!/bin/bash
#SBATCH --job-name=stageEad
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr4
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --array=0-4
#SBATCH --output=/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/stage_e/stage_e_adaptive_%a_%j.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/stage_e/stage_e_adaptive_%a_%j.err
# One alpha per array task (0..4 -> alpha 0.40,0.45,0.50,0.55,0.60), parallel.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
GWENV=/global/scratch/projects/pc_heptheory/clwelch/software/envs/gwenv/bin/python
cd /global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d/stage_e
$GWENV stage_e_adaptive.py
