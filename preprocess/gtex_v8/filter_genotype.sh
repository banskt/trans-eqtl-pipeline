#!/bin/bash -e

INVCF="/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_chr[CHRM].vcf.gz"
OUTVCF="/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_chr[CHRM].vcf.gz"

# Define all preprocessing scripts
PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
PYFILTER="${PWD}/scripts/filter_vcf.py"

source "${PWD}/utils/submit_job"

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PWD}/jobsubs"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
JOBNAME="gtex_v8_gtfilter_[CHRM]_${RANDSTRING}"

for CHRM in `seq 1 22`; do 
    CHRM_INVCF=${INVCF/\[CHRM\]/${CHRM}}
    CHRM_OUTVCF=${OUTVCF/\[CHRM\]/${CHRM}}
    CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

    sed -e "s|_PYENV_|${PYENV}|g;
            s|_PYFILTER_|${PYFILTER}|g;
            s|_JOB_NAME|${CHRM_JOBNAME}|g;
            s|_INVCF_|${CHRM_INVCF}|g;
            s|_OUTVCF_|${CHRM_OUTVCF}|g
            " ${PWD}/bsubfiles/vcf_filter.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

    submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
done