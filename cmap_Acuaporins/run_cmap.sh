#!/bin/bash

# limitar threads internos
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# opcional: controlar cores del script
#export N_CORES=6

# ejecutar script
nohup python3.9 cmap.py &
