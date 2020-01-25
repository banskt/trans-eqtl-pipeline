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

PYPHASER="${PREPROC_SCRIPTDIR}/process_phASER_counts.py"
PYTPMS="${PREPROC_SCRIPTDIR}/calculate_TPMs.py"

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

if [ "${USE_BSUB}" = "true" ]; then
    source utils/preprocess_reads_and_tpms.sh
else
    if [ "${bProcessPhaser}" = "true" ]; then source ${PREPROC_UTILSDIR}/000_process_phASER_counts.sh; fi
    if [ "${bCalculateTPMs}" = "true" ]; then source ${PREPROC_UTILSDIR}/001_process_TPMs_from_reads.sh; fi
fi