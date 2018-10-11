#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./01_validation_pipeline.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

for MDATA in ${DATASETS}; do

    source PATHS
    JOBSUBDIR_DATA="${JOBSUBDIR}/${MDATA}"
    OUTDIR_DATA="${OUTDIR}/${MDATA}"

    if [ ! -z "$EXPRESSIONFILE" ]; then

        if [ "${bMatrixEqtl}" = "true" ]; then source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bTejaas}" = "true" ]; then source ${UTILSDIR}/tejaas; fi
        echo "Submitting jobs for $MDATA"

    fi

    source unset_variables.sh PATHS

done

if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

source unset_variables.sh ${CONFIGFILE}
