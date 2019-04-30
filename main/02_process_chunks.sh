#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_process_chunks.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/gx_preproc_string
source ${UTILSDIR}/tejaas_chunk_reduce.new

source ${DATALOAD}

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        OUTDIR_DATA="${OUTDIR}/${TSHORT}"
        EXPRESSIONFILE=${EXPR_FMT/\[TISSUE\]/${TSHORT}}

        echo ${OUTDIR_DATA}

        if [ ! -z "$EXPRESSIONFILE" ]; then

            echo "Processing chunks for $TSHORT"

            for CHRM in ${CHRNUMS}; do
                # calculate total number of chunks expected
                NCHUNK=$( expected_nchunk ${GENO_FMT} ${CHRM} ${MAX_NSNP_PERJOB} )
                echo "chr${CHRM} --> ${NCHUNK} chunks"
                echo "============================"
                if [ "${bTejaasJPA}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/jpa/chr${CHRM}" $NCHUNK; fi
                if [ "${bJPARandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/jpa/chr${CHRM}" $NCHUNK; fi
                for NULL in ${TEJAAS_NULL}; do
                    if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                    if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                    for SBETA in ${TEJAAS_SIGMA_BETA}; do
                        if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${NULL}null_sb${SBETA}/chr${CHRM}" $NCHUNK; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${NULL}null_sb${SBETA}/chr${CHRM}" $NCHUNK; fi
                    done
                done
            done
        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
