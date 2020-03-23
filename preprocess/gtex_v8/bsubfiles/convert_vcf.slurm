#!/bin/sh
#SBATCH -p hh
#SBATCH -t 1-00:00:00
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174

PYENV=_PYENV_
INVCF=_INVCF_
OUTDS=_OUTDS_
PYCONVERT=_PYCONVERT_

$PYENV $PYCONVERT --in $INVCF \
                 --out $OUTDS 