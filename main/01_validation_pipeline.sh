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
    SHUFFLED_ID_FILE="${OUTDIR}/shuffled_donor_ids.txt"

    source ${UTILSDIR}/shuffle_donors
    

    if [ ! -z "$EXPRESSIONFILE" ]; then

        echo "Submitting jobs for $MDATA"
        
        if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bMEqtlRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; fi
        if [ "${bTejaas}" = "true" ];      then source ${UTILSDIR}/tejaas; fi
        if [ "${bTjsRandom}" = "true" ];   then SHUFFLE=true; source ${UTILSDIR}/tejaas; fi
        if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi
        if [ "${bJPARandom}" = "true" ];   then SHUFFLE=true; RUNJPA=true; source ${UTILSDIR}/tejaas; fi

        if [ "${bTjsRandom1000}" = "true" ];   then SHUFFLE=true; source ${UTILSDIR}/tejaas_randoms; fi

        if [ "${bTejaasPartition}" = "true" ]; then
            for r in `seq 1 $REPLICAS`; do
                echo $r
                PARTITION_FILE_1="${OUTDIR_DATA}/${r}/part1.txt"
                PARTITION_FILE_2="${OUTDIR_DATA}/${r}/part2.txt"
                
                # Create partitioned donor files
                source ${UTILSDIR}/partition_donors

                PARTITION1=true; PARTITION2=false;
                source ${UTILSDIR}/tejaas;
                # source ${UTILSDIR}/matrix_eqtl;
                # SHUFFLE=true; source ${UTILSDIR}/tejaas
                # SHUFFLE=false

                PARTITION1=false; PARTITION2=true;
                source ${UTILSDIR}/tejaas;
                # source ${UTILSDIR}/matrix_eqtl;
                # SHUFFLE=true; source ${UTILSDIR}/tejaas
                # SHUFFLE=false
            done;
        fi
    fi

    # echo ${JOBDEPS}

    unset_vars DATA

done

if [ "${bValidationPlot}" = "true" ]; then source ${UTILSDIR}/validation_plot; fi

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
