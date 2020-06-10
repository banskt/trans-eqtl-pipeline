#!/bin/bash
#SBATCH -p hh
#SBATCH -t 2-0:00:00
#SBATCH --exclusive
#SBATCH -J sboptim_fram_all_30
#SBATCH -o sboptim_fram_all_30.out
#SBATCH -e sboptim_fram_all_30.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
export OMP_NUM_THREADS=16
MPYTHON=/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python

RUN_PATH=/usr/users/fsimone/saikat/sb_optimization/sigbeta_optimizer.py
GXPRFILE=/cbscratch/franco/datasets/FHS/expression/fhs.formatted.expr.one_pedno_rep.txt
OUTFILE="/usr/users/fsimone/saikat/sb_optimization/framingham_uniq_sbpath_30.txt"

${MPYTHON} ${RUN_PATH} --gx ${GXPRFILE} --out ${OUTFILE} --nsnp 5000
