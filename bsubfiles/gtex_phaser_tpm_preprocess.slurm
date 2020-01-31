#!/bin/sh
#SBATCH -p hh
#SBATCH -t 0-2:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
# module load R/3.4.1
PYENV=_PYENV_
PYPHASER=_PYPHASER_
PYTPMS=_PYTPMS_
PREPROC_UTILSDIR=_PUTILS_DIR_

# export OMP_NUM_THREADS=8

bUsePub=_USE_PUB_

# source files
GTFFILE=_GTFFILE_
OUTTPM=_OUTTPM_
SRCREAD=_SRCREAD_
PHASER_INMTX=_PHASER_IN_
PHASER_OUTMTX=_PHASER_OUT_

# which steps to be run
bProcessPhaser=_STEP1_
bCalculateTPMs=_STEP2_

# Select tissue-specific rpkms and normalize
echo "01 - Process phASER counts"
if [ "${bProcessPhaser}" = "true" ]; then source ${PREPROC_UTILSDIR}/000_process_phASER_counts.sh; fi

# Obtain Covariates from GTEX
echo "02 - Calculate TPMs from counts"
if [ "${bCalculateTPMs}" = "true" ]; then source ${PREPROC_UTILSDIR}/001_process_TPMs_from_reads.sh; fi
