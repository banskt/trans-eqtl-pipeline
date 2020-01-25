#!/bin/bash -e
source ${UTILSDIR}/submit_job

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PREPROCDIR}/gtex_v8/jobsubs"
# if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

JOBNAME="gtex_phaser_tpms_preprocess_gx_${RANDSTRING}"

sed -e "s|_JOB_NAME|${JOBNAME}|g;
        s|_PYENV_|${PYENV}|g;
        s|_PYPHASER_|${PYPHASER}|g;
        s|_PYTPMS_|${PYTPMS}|g;
        s|_OUTTPM_|${OUTTPM}|g;
        s|_SRCREAD_|${SRCREAD}|g;
        s|_GTFFILE_|${GTFFILE}|g;
        s|_STEP1_|\"${bProcessPhaser}\"|g;
        s|_STEP2_|\"${bCalculateTPMs}\"|g;
        s|_USE_PUB_|${bUsePub}|g;
        s|_PHASER_IN_|${PHASER_INMTX}|g;
        s|_PHASER_OUT_|${PHASER_OUTMTX}|g;
        s|_PUTILS_DIR_|${PREPROC_UTILSDIR}|g;
        " ${MASTER_BSUBDIR}/gtex_phaser_tpm_preprocess.slurm > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.slurm

submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
   

