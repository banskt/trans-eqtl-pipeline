#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./01_validation_pipeline.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_chunk_reduce

for MDATA in ${DATASETS}; do

    source DATA
    OUTDIR_DATA="${OUTDIR}/${MDATA}"

    if [ ! -z "$EXPRESSIONFILE" ]; then

        echo "Processing chunks for $MDATA"

        for CHRM in ${CHRNUMS}; do
            if [ "${bTejaasJPA}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/jpa/chr${CHRM}"; fi
            if [ "${bJPARandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/jpa/chr${CHRM}"; fi
            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for SBETA in ${TEJAAS_SIGMA_BETA}; do
                    if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${NULL}null_sb${SBETA}/chr${CHRM}"; fi
                    if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${NULL}null_sb${SBETA}/chr${CHRM}"; fi
                done
            done
        done
    fi

    unset_vars DATA

done

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
