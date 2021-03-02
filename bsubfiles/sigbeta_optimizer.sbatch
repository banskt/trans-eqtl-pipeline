#!/bin/bash
#SBATCH -A cramer
#SBATCH -p em
#SBATCH -t 2-0:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
# module load openmpi/intel/64/1.10.7

export OMP_NUM_THREADS=8
MPYTHON=_PYT_ENV_

RUN_PATH=_SBOP_PY_
GXPRFILE=_EXPR_FL_
SAMPFILE=_SAMP_FL_
OUTFILE=_OUTF_SB_

${MPYTHON} ${RUN_PATH} --gx ${GXPRFILE} --fam ${SAMPFILE} --out ${OUTFILE} --nsnp 5000 --knn 30
