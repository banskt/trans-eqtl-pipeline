#!/bin/bash
source ${UTILSDIR}/submit_job

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PREPROCDIR}/framingham/jobsubs"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

for CONSENT in ${CONSENT_GROUPS}; do
    SRCDOSE_CNST=${SRCDOSE//\[CNST\]/${CONSENT}}
    GTOUT_CNST=${GT_OUTFILE//\[CNST\]/${CONSENT}}

    for CHRM in ${CHRNUMS}; do
        JOBNAME="fhs_gt_preprocess_${CONSENT}_${CHRM}_${RANDSTRING}"
        GENOTYPEFILE=${SRCDOSE_CNST/\[CHRM\]/${CHRM}}
        SRCINFO_CHRM=${SRCINFO/\[CHRM\]/${CHRM}}
        DOSE_OUTFILE=${GTOUT_CNST/\[CHRM\]/${CHRM}}

        if [ -s ${DOSE_OUTFILE} ]; then
            echo "File exists: ${DOSE_OUTFILE}"
            exit;
        fi
        
        sed -e "s|_JOB_NAME|${JOBNAME}|g;
                s|_RFSCRIPT_|${REFORMATPY}|g;
                s|_PYENV_|${PYENV}|g;
                s|_GTOUT_|${DOSE_OUTFILE}|g;
                s|_SMPLE_|${SRCSAMPLE}|g;
                s|_CNST_|${CONSENT}|g;
                s|_GTIN_|${GENOTYPEFILE}|g;
                s|_INFO_|${SRCINFO_CHRM}|g;
                s|_CHRM_|${CHRM}|g;
                s|_KEEP_|${STARTERSNPS}|g;
                " ${MASTER_BSUBDIR}/fhs_gt_preprocess.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

                # s|_GX_PROC_SH_|${UTILSDIR}/gx_preproc_string|g;
                # s|_GX_FIL_FMT_|${GXFILENAME_FMT}|g;

        submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}

    done;
done;

