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

#---- Include functions
# source ${UTILSDIR}/gx_preproc_string
# source ${UTILSDIR}/submit_job

# Define all preprocessing scripts
PREPROC_SCRIPTDIR="${PWD}/scripts"
PREPROC_UTILSDIR="${PWD}/utils"
REFORMATPY="${PREPROC_SCRIPTDIR}/reformat_dosages.py"
CONVERTPY="${PREPROC_SCRIPTDIR}/dosages_to_vcf.py"

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

if [ "${bReformatGT}" = "true" ]; then source ${PREPROC_UTILSDIR}/reformat_gt.sh; fi
if [ "${bMergeCnst}" = "true" ]; then source ${PREPROC_UTILSDIR}/merge_gt_groups.sh; fi
if [ "${bConvert2VCF}" = "true" ]; then source ${PREPROC_UTILSDIR}/convert2vcf.sh; fi