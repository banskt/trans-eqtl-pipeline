#!/bin/bash

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/preprocess/gtex/expression"
if [ ! -d ${SPECIFIC_JOBSUBDIR} ]; then mkdir -p ${SPECIFIC_JOBSUBDIR}; fi

GX_TISSUE_JOBNAME="gtex_preprocess_gx_${TSHORT}_${RANDSTRING}"

if [ "${USE_SLURM}" = "true" ]; then
    __MASTERFILE="slurmfiles/gx_preproc.sbatch"
    __SUBMITFILE="${SPECIFIC_JOBSUBDIR}/${GX_TISSUE_JOBNAME}.sbatch"
elif [ "${USE_LSF}" = "true" ]; then
    __MASTERFILE="lsffiles/gx_preproc.bsub"
    __SUBMITFILE="${SPECIFIC_JOBSUBDIR}/${GX_TISSUE_JOBNAME}.bsub"
fi

sed -e "s|_PYENV_|${PYENV}|g;
        s|_JOB_NAME|${GX_TISSUE_JOBNAME}|g;
        s|_T_FULL_|\"${TFULL}\"|g;
        s|_T_SHORT_|${TSHORT}|g;
        s|_T_BASE_|${TBASE}|g;
        s|_SRC_TPM_|${SRCTPM}|g;
        s|_SRCREAD_|${SRCREAD}|g;
        s|_SRCPHENO_|${SRCPHENO}|g;
        s|_SRCCOVARS_|${SRCCOVARS}|g;
        s|_SRCGENCODE_|${GENCODEFILE}|g;
        s|_PUTILS_DIR_|${PREPROC_UTILSDIR}|g;
        s|_OUT_DIR_|${OUTDIR}|g;
        s|_TIS_OUT_|${TISSUEOUTDIR}|g;
        s|_PRE_GXO_|${PREGXOUTDIR}|g;
        s|_GXF_OUT_|${GXOUTDIR}|g;
        s|_COV_DIR_|${COVARDIR}|g;
        s|_NORM_QC_PY_|${NORMQCPY}|g;
        s|_TPM_THR_|${TPM_THRESHOLD}|g;
        s|_CNT_THR_|${COUNTS_THRESHOLD}|g;
        s|_SAM_FRX_|${SAMPLE_FRAC_THRESHOLD}|g;
        s|_QC_MTHD_|\"${QCMETHODS}\"|g;
        s|_GXT_SEL_|\"${GXSELECTION}\"|g;
        s|_b_CLCOV_|${bCollectCovs}|g;
        s|_b_SELTS_|${bSelectTissue}|g;
        s|_b_NORQC_|${bNormalizeQC}|g;
    " ${__MASTERFILE} > ${__SUBMITFILE}

if [ "${USE_SLURM}" = "true" ]; then
    GX_TISSUE_JOBID=$( submit_job ${SPECIFIC_JOBSUBDIR} ${GX_TISSUE_JOBNAME} ${THISJOBDEPS} )
    GX_TISSUE_JOBDEPS=$( add_deps "${GX_TISSUE_JOBDEPS}" ${GX_TISSUE_JOBID} )
elif [ "${USE_LSF}" = "true" ]; then
    GX_TISSUE_JOBID=$( submit_job_lsf ${SPECIFIC_JOBSUBDIR} ${GX_TISSUE_JOBNAME} ${THISJOBDEPS} )
    GX_TISSUE_JOBDEPS=$( add_deps_lsf "${GX_TISSUE_JOBDEPS}" ${GX_TISSUE_JOBID} )
fi
