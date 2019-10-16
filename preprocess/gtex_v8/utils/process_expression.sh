#!/bin/bash -e
source ${UTILSDIR}/submit_job

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PREPROCDIR}/gtex_v8/jobsubs"
# if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

JOBNAME="gtex_preprocess_gx_${TSHORT}_${RANDSTRING}"

sed -e "s|_JOB_NAME|${JOBNAME}|g;
        s|_PYENV_|${PYENV}|g;
        s|_T_FULL_|\"${TFULL}\"|g;
        s|_T_SHORT_|\"${TSHORT}\"|g;
        s|_T_BASE_|\"${TBASE}\"|g;
        s|_SRCTPM_|${SRCTPM}|g;
        s|_SRCPHENO_|${SRCPHENO}|g;
        s|_SRCREAD_|${SRCREAD}|g;
        s|_SRCCOVARS_|${SRCCOVARS}|g;
        s|_SRCGENCODE_|${GENCODEFILE}|g;
        s|_DONOR_F_|${DONORFILE}|g;
        s|_GX_OUT_|${GXOUTDIR}|g;
        s|_LM_OUT_|${LMOUTDIR}|g;
        s|_PEER_OUT_|${PEEROUTDIR}|g;
        s|_COV_OUT_|${COVOUTDIR}|g;
        s|_TPM_OUT_|${TPMOUTDIR}|g;
        s|_AGECOV_F_|${AGE_COVARIATE_FILE}|g;
        s|_CORRAGE_|\"${GXAGECORR}\"|g;
        s|_STEP1_|\"${bSelectNormalize}\"|g;
        s|_STEP2_|\"${bFormatCovariates}\"|g;
        s|_STEP3_|\"${bLMcorrect}\"|g;
        s|_STEP4_|\"${bPeerCorrect}\"|g;
        s|_STEP5_|\"${bGencodeFilter}\"|g;
        s|_NCOVS_|\"${GXNPEERS}\"|g;
        s|_GX_SEL_|\"${GXSELECTION}\"|g;
        s|_SL_SAMP_PY_|${SELECTSAMPLEPY}|g;
        s|_GX_NORM_PY_|${GTEXNORMALIZEPY}|g;
        s|_AGE_COV_PY_|${COMPILEAGECOVPY}|g;
        s|_LM_CORR_RS_|${LMCORR_R}|g;
        s|_CORR_PY_|${CORRECTPY}|g;
        s|_PR_CORR_RS_|${PEERSCRIPT_R}|g;
        s|_GEN_FIL_PY_|${GENCODEFILTERPY}|g;
        s|_PUTILS_DIR_|${PREPROC_UTILSDIR}|g;
        s|_USE_PUB_|${bUsePub}|g;
        " ${MASTER_BSUBDIR}/gtex_gx_preprocess.slurm > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.slurm

        # " ${MASTER_BSUBDIR}/gtex_gx_preprocess.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

        # s|_GX_PROC_SH_|${UTILSDIR}/gx_preproc_string|g;
        # s|_GX_FIL_FMT_|${GXFILENAME_FMT}|g;

submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
# submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
