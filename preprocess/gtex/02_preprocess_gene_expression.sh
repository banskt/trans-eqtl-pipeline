#!/bin/bash -e

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi
CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi
source ${CONFIGFILE}
source ${DATALOAD}
source ${EXTERNALLOAD}
source ../../main/PATHS

#---- Include functions from main directory
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps

PREPROC_SCRIPTDIR="${PWD}/scripts"
PREPROC_UTILSDIR="${PWD}/utils"
NORMQCPY="${PREPROC_SCRIPTDIR}/preprocess_expression.py"

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

# Create list of samples
zcat ${SRCTPM} | head -n 3 | tail -n 1 | sed 's/\t/\n/g' > ${PREGXOUTDIR}/rnaseq_samples.txt
cut -d" " -f1 ${DONORFILE} | tail -n +3 > ${PREGXOUTDIR}/vcf_samples.list

#Collect covariates
source ${PREPROC_UTILSDIR}/collect_age_gender_trischd.sh

## control job dependencies
GX_TISSUE_JOBDEPS="None"

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        ### Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
        TISSUEOUTDIR="${PREGXOUTDIR}/${TSHORT}"
        if [ ! -d ${TISSUEOUTDIR} ]; then mkdir -p ${TISSUEOUTDIR}; fi

        echo $TFULL

        if [ "${USE_SLURM}" = "true" ]; then
            echo "Submitting jobs using SLURM"
            source ${PREPROC_UTILSDIR}/submit_gx_preproc.sh
        elif [ "${USE_LSF}" = "true" ]; then
            echo "Submitting jobs using LSF is no longer supported. Please move to SLURM"
            ## source ${PREPROC_UTILSDIR}/submit_process_expression.sh
        else
            echo "Doing it serially."
            source ${PREPROC_UTILSDIR}/gx_preproc.sh
        fi
    fi
done < ${TISSUEFILE}

echo ${GX_TISSUE_JOBDEPS}
