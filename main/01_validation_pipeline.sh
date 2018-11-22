#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./01_validation_pipeline.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/unset_vars
source EXTERNAL

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
RUNJPA=false	# used for submitting jpa-only jobs
SHUFFLE=false 	# used for controlling shuffling
JOBDEPS="None" 	# used for controlling job dependencies

for MDATA in ${DATASETS}; do

    source DATA
    JOBSUBDIR_DATA="${JOBSUBDIR}/${MDATA}"
    OUTDIR_DATA="${OUTDIR}/${MDATA}"

    if [ ! -z "$EXPRESSIONFILE" ]; then

        echo "Submitting jobs for $MDATA"
        
        if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bMEqtlRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bTejaas}" = "true" ];      then source ${UTILSDIR}/tejaas; fi
        if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi

    fi

    # echo ${JOBDEPS}

    unset_vars DATA

done

if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

unset_vars ${CONFIGFILE}
unset_vars PATHS
