#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./01_validation_pipeline.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

# get the different functions
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/unset_vars

for MDATA in ${DATASETS}; do

    source DATA
    JOBSUBDIR_DATA="${JOBSUBDIR}/${MDATA}"
    OUTDIR_DATA="${OUTDIR}/${MDATA}"

    if [ ! -z "$EXPRESSIONFILE" ]; then

        echo "Submitting jobs for $MDATA"
        JOBDEPS=""
        
        if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bMEqtlRandom}" = "true" ]; then source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bTejaas}" = "true" ];      then source ${UTILSDIR}/tejaas; fi

    fi

    unset_vars DATA

done

if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

unset_vars ${CONFIGFILE}
unset_vars PATHS
