#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_process_chunks.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_chunk_reduce.new

source ${DATALOAD}

## chromosomes in separate outer loop to prevent reading VCF file multiple times.
NCHUNK_IN_CHRM=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # initialize with 22 elements
for CHRM in ${CHRNUMS}; do 
    CHRINDX=$((CHRM - 1))
    SPECIFIC_OUTDIR="/usr/users/sbanerj/trans-eQTL/dev-pipeline/jobsubs/gtex_v8/as/tejaas/raw_std/permnull_sb0.1_knn/chr${CHRM}"
    NCHUNK=$( ls -l ${SPECIFIC_OUTDIR}/*.sbatch | wc -l )
    #NCHUNK=$( expected_nchunk ${GENO_FMT} ${CHRM} ${MAX_NSNP_PERJOB} )
    NCHUNK_IN_CHRM[${CHRINDX}]=${NCHUNK}
    echo "Chr${CHRM} was split into ${NCHUNK} jobs."
done

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        OUTDIR_DATA="${OUTDIR}/${TSHORT}"
        GX_TISSUE_FMT=${EXPR_FMT/\[TISSUE\]/${TSHORT}}
        EXPRESSIONFILE=${GX_TISSUE_FMT/\[PREPROC_STRING\]/${TEJAAS_PREPROC_STR}}

        echo ${OUTDIR_DATA}

        if [ ! -z "$EXPRESSIONFILE" ]; then

            echo "Processing chunks for $TSHORT"

            for CHRM in ${CHRNUMS}; do
                CHRINDX=$((CHRM - 1))
                NCHUNK=${NCHUNK_IN_CHRM[${CHRINDX}]}
                echo "chr${CHRM} --> ${NCHUNK} chunks"
                echo "============================"
                for NULL in ${TEJAAS_NULL}; do
                    if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                    if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                    for SBETA in ${TEJAAS_SIGMA_BETA}; do
                        METHOD_VARIANT="${NULL}null_sb${SBETA}"
                        if [ "${TEJAAS_KNN}" = "true" ]; then METHOD_VARIANT="${METHOD_VARIANT}_knn"; fi
                        if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}" $NCHUNK; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}" $NCHUNK; fi
                    done
                done
            done
        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
