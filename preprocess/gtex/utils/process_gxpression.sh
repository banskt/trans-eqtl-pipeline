#!/bin/bash

THISMETHOD="process_gx"
JOBNAME="${THISMETHOD}_${TSHORT}"

SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/${THISMETHOD}"

# if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR};

sed -e "s|_T_FULL_|\"${TFULL}\"|g;
        s|_T_SHORT_|\"${TSHORT}\"|g;
        s|_T_BASE|\"${TBASE}\"|g;
        s|_UTILS_|${UTILSDIR}|g;
        s|_PYENV_|${PYENV}|g;
        s|_SCRIPTDIR_|${SCRIPTDIR}|g;
        s|_SRCRPKM_|${SRCRPKM}|g;
        s|_SRCPHENO_|${SRCPHENO}|g;
        s|_SRCREAD_|${SRCREAD}|g;
        s|_SRCCOVARS_|${SRCCOVARS}|g;
        s|_DONORF_|${DONORFILE}|g;
        s|_NORM_OUTDIR_|${NORMOUTDIR}|g;
        s|_LM_OUT_|${LMOUTDIR}|g;
        s|_COVOUTDIR_|${COVOUTDIR}|g;
        s|_RPKM_OUT_|${RPKMOUTDIR}|g;
        s|_NCOVS_|\"${NCOVS}\"|g;
        s|_STEP1_|\"${bProcessRPKMs_and_normalize}\"|g;
        s|_STEP2_|\"${bformatCovariates}\"|g;
        s|_STEP3_|\"${bLMcorrect}\"|g;
        s|_STEP4_|\"${GTEX_PEER_CORRECTION}\"|g;
        s|_CORRAGE_|\"${CORRECT_AGE}\"|g;
        s|_AGE_COV_|${AGE_COVARIATE}|g;
        s|_PEEROUT_|${PEEROUTDIR}|g;
        s|_JOB_NAME|${JOBNAME}|g;
        " ${MASTER_BSUBDIR}/process_gx.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} None # no job dependencies

