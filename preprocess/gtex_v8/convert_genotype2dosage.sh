#!/bin/bash -e

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi
CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source "${PWD}/utils/submit_job"

for MAF_THRES in $MAF_THRESHOLDS; do
    THIS_VCFDIR="${INVCFDIR}/${MAF_THRES}"
    INVCF="${THIS_VCFDIR}/${VCFOBASE}.vcf.gz"
    OUTDS="${THIS_VCFDIR}/${VCFOBASE}.dosage.gz"



    THISJOBDEPS="None" # no job dependencies
    SPECIFIC_JOBSUBDIR="${PWD}/jobsubs/genotype"
    #if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

    # Define a random string for marking jobs of this batch
    RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
    JOBNAME="gtex_v8_vcf2dosage_[CHRM]_${RANDSTRING}"

    for CHRM in `seq 1 22`; do 
        CHRM_INVCF=${INVCF/\[CHRM\]/${CHRM}}
        CHRM_OUTDS=${OUTDS/\[CHRM\]/${CHRM}}
        CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

        CHRM_INVCF=${CHRM_INVCF/\[MAFTHRES\]/${MAF_THRES}}
        CHRM_OUTDS=${CHRM_OUTDS/\[MAFTHRES\]/${MAF_THRES}}

        sed -e "s|_PYENV_|${PYENV}|g;
                s|_PYCONVERT_|${PYCONVERT}|g;
                s|_JOB_NAME|${CHRM_JOBNAME}|g;
                s|_INVCF_|${CHRM_INVCF}|g;
                s|_OUTDS_|${CHRM_OUTDS}|g
                " ${PWD}/bsubfiles/convert_vcf.slurm > ${SPECIFIC_JOBSUBDIR}/${CHRM_JOBNAME}.slurm

        submit_job ${SPECIFIC_JOBSUBDIR} ${CHRM_JOBNAME} ${THISJOBDEPS}
    done
done