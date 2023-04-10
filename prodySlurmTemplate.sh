#!/bin/bash
#SBATCH	-n 6
#SBATCH	--job-name=FaPIP21_dih
#SBATCH	--partition=intel
#SBATCH	--no-requeue
#SBATCH	--mem-per-cpu=1G

export OMP_NUM_THREADS=$SLURM_NPROCS
python script.py [extra_arguments]
