#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/gx_preproc_string
source EXTERNAL

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
RUNJPA=false	# used for submitting jpa-only jobs
SHUFFLE=false 	# used for controlling shuffling
JOBDEPS="None" 	# used for controlling job dependencies
TEJAAS_JOBDEPS="None"
SUBMITTED_JOBIDS="" # used for controlling jobid reporting

source ${DATALOAD}

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        JOBSUBDIR_DATA="${JOBSUBDIR}/${TSHORT}"
        OUTDIR_DATA="${OUTDIR}/${TSHORT}"
        SHUFFLED_ID_FILE="${OUTDIR}/shuffled_donor_ids.txt"

        GX_TISSUE_FMT=${EXPR_FMT/\[TISSUE\]/${TSHORT}}
        EXPRESSIONFILE=${GX_TISSUE_FMT/\[PREPROC_STRING\]/${TEJAAS_PREPROC_STR}}

        source ${UTILSDIR}/shuffle_donors

        if [ ! -z "$EXPRESSIONFILE" ]; then

            echo "Submitting jobs for $TBASE"
        
            if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
            if [ "${bMEqtlRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; SHUFFLE=false; fi
            if [ "${bTejaas}" = "true" ];      then source ${UTILSDIR}/tejaas; fi
            if [ "${bTjsRandom}" = "true" ];   then SHUFFLE=true; source ${UTILSDIR}/tejaas; SHUFFLE=false; fi
            #if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi
            #if [ "${bJPARandom}" = "true" ];   then SHUFFLE=true; RUNJPA=true; source ${UTILSDIR}/tejaas; fi
            #if [ "${bGNetLmm}" = "true" ]; then source ${UTILSDIR}/gnetlmm; fi

        fi
    # echo ${JOBDEPS}
    fi
done < ${TISSUEFILE}

#if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
