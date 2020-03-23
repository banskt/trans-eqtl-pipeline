#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/unset_vars
#source ${UTILSDIR}/gx_preproc_string
PYMQTL="${SCRIPTDIR}/reannotate_meqtl_results.py"

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
RUNJPA=false	# used for submitting jpa-only jobs
SHUFFLE=false 	# used for controlling shuffling
JOBDEPS="None" 	# used for controlling job dependencies

for EXPR_CORR in ${EXPRESSIONS}; do
    for MDATA in ${DATASETS}; do
        # source DATA
        source ${DATALOAD}
        JOBSUBDIR_DATA="${JOBSUBDIR}/${MDATA}/$EXPR_CORR"
        echo $JOBSUBDIR_DATA
        OUTDIR_DATA="${OUTDIR}/${MDATA}"
        if [ ! -d ${OUTDIR_DATA} ]; then mkdir -p ${OUTDIR_DATA}; fi
        SHUFFLED_ID_FILE="${OUTDIR_DATA}/shuffled_donor_ids.txt"

        source ${UTILSDIR}/shuffle_donors

        if [ ! -z "$EXPRESSIONFILE" ]; then

            echo "Submitting jobs for $MDATA"
            
            if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
            if [ "${bMEqtlRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; fi

            for KNN in $KNNS; do
                if [ "${bTejaas}" = "true" ];      then source ${UTILSDIR}/tejaas; fi
                if [ "${bTjsRandom}" = "true" ];   then SHUFFLE=true; source ${UTILSDIR}/tejaas; fi
                if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi
                if [ "${bJPARandom}" = "true" ];   then SHUFFLE=true; RUNJPA=true; source ${UTILSDIR}/tejaas; fi


                if [ "${bTjsRandomN}" = "true" ];   then SHUFFLE=true; source ${UTILSDIR}/tejaas_randoms; fi
            done
        fi
        
        # echo ${JOBDEPS}
        unset_vars ${DATALOAD}

    done
done

# if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

unset_vars PATHS
unset_vars ${EXTERNALLOAD}
unset_vars ${CONFIGFILE}