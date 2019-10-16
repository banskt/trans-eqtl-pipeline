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

# export OMP_NUM_THREADS=8

# tissue name
TFULL=_T_FULL_
TSHORT=_T_SHORT_
TBASE=_T_BASE_
bUsePub=_USE_PUB_

# source files
SRCRPKM=_SRCRPKM_
SRCTPM=_SRCTPM_
SRCPHENO=_SRCPHENO_
SRCREAD=_SRCREAD_
SRCCOVARS=_SRCCOVARS_
GENCODEFILE=_SRCGENCODE_
DONORFILE=_DONOR_F_

# output directories / files
GXOUTDIR=_GX_OUT_
LMOUTDIR=_LM_OUT_
PEEROUTDIR=_PEER_OUT_
COVOUTDIR=_COV_OUT_
RPKMOUTDIR=_RPKM_OUT_
TPMOUTDIR=_TPM_OUT_
AGE_COVARIATE_FILE=_AGECOV_F_

# which steps to be run
GXAGECORR=_CORRAGE_
bSelectNormalize=_STEP1_
bFormatCovariates=_STEP2_
bLMcorrect=_STEP3_
bPeerCorrect=_STEP4_
bGencodeFilter=_STEP5_

# parameters for preprocessing
GXNPEERS=_NCOVS_
GXSELECTION=_GX_SEL_

# pipeline scripts / directories
SELECTSAMPLEPY=_SL_SAMP_PY_
GTEXNORMALIZEPY=_GX_NORM_PY_
COMPILEAGECOVPY=_AGE_COV_PY_
CORRECTPY=_CORR_PY_
# LMCORR_R=_LM_CORR_RS_
PEERSCRIPT_R=_PR_CORR_RS_
GENCODEFILTERPY=_GEN_FIL_PY_
# GX_PREPROC_STR_BASH=_GX_PROC_SH_
PREPROC_UTILSDIR=_PUTILS_DIR_

# final output filename format
# GXFILENAME_FMT=_GX_FIL_FMT_

# source ${GX_PREPROC_STR_BASH}

# Expression Files
# NORMOUTFILE="${NORMOUTDIR}/${TSHORT}_normalized.txt"
# LMOUTFILE="${LMOUTDIR}/${TSHORT}_lmcorrected.txt"
# LMOUTFILE_AGE="${LMOUTDIR}/${TSHORT}_age_lmcorrected.txt"

# Covariate Files
COVARS="${COVOUTDIR}/${TSHORT}_nopeer_covariates.txt"
COVARS_AGE="${COVOUTDIR}/${TSHORT}_nopeer_covariates_w_age.txt"

# Select tissue-specific rpkms and normalize
echo "01 - QC and Normalization"
if [ "${bSelectNormalize}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_01_select_normalize.sh; fi

# Obtain Covariates from GTEX
echo "02 - Formating Covariates"
if [ "${bFormatCovariates}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_02_format_covariates.sh; fi

# Correct gene expression by linear model regression of covariates
echo "03 - Linear Model Correction"
if [ "${bLMcorrect}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_03_lmcorrect.sh; fi 

# Correct gene expression by PEER factors
echo "04 - PEER Correction"
if [ "${bPeerCorrect}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_04_peer.sh; fi 
 
# Filter with GENCODE
echo "05 - GENCODE filtering"
if [ "${bGencodeFilter}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_05_gencode_filter.sh; fi
