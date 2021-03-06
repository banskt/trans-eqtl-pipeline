#!/bin/bash
#SBATCH -p medium
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -C scratch
#SBATCH -A all
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load R/3.4.1

export OMP_NUM_THREADS=2

## For required variables, see the file ${PREPROC_UTILSDIR}/gx_preproc.sh

PYENV=_PYENV_

# tissue name
TFULL=_T_FULL_
TSHORT=_T_SHORT_
TBASE=_T_BASE_

# source files
SRCTPM=_SRC_TPM_
SRCREAD=_SRCREAD_
SRCPHENO=_SRCPHENO_
SRCCOVARS=_SRCCOVARS_
GENCODEFILE=_SRCGENCODE_

# pipeline scripts / directories
PREPROC_UTILSDIR=_PUTILS_DIR_
OUTDIR=_OUT_DIR_
GXOUTDIR=_GXF_OUT_
PREGXOUTDIR=_PRE_GXO_
TISSUEOUTDIR=_TIS_OUT_
COVARDIR=_COV_DIR_
NORMQCPY=_NORM_QC_PY_

# parameters for preprocessing
TPM_THRESHOLD=_TPM_THR_
COUNTS_THRESHOLD=_CNT_THR_
SAMPLE_FRAC_THRESHOLD=_SAM_FRX_
QCMETHODS=_QC_MTHD_
GXSELECTION=_GXT_SEL_

bCollectCovs=_b_CLCOV_
bSelectTissue=_b_SELTS_
bNormalizeQC=_b_NORQC_

source ${PREPROC_UTILSDIR}/gx_preproc.sh
