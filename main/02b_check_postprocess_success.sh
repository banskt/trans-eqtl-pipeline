#!/bin/bash
if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars

source ${DATALOAD}

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/postprocess"
        OUTDIR_DATA="${OUTDIR}/${TSHORT}"

        GX_TISSUE_FMT=${EXPR_FMT/\[TISSUE\]/${TSHORT}}
        EXPRESSIONFILE=${GX_TISSUE_FMT/\[PREPROC_STRING\]/${TEJAAS_PREPROC_STR}}

        if grep -Fwq ${TSHORT} ${TEJAAS_SIGMA_BETA_PERM_FILE}; then
            TEJAAS_SIGMA_BETA_PERM=$( grep -w ${TSHORT} ${TEJAAS_SIGMA_BETA_PERM_FILE} | cut -f3 )
        fi

        if [ ! -z "$EXPRESSIONFILE" ]; then

            METHOD_VARIANT="${TEJAAS_NULL}null_sb${TEJAAS_SIGMA_BETA_PERM}"
            if [ "${TEJAAS_KNN}" = "true" ]; then
                METHOD_VARIANT="${METHOD_VARIANT}_knn${KNN_NBR}"
            fi

            RESULT_DIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}"
            COLRES_DIR="${TEQTL_OUTDIR}/${TSHORT}"

            failedChunkReduce=false
            for CHRM in ${CHRNUMS}; do
                for EXT in rr.txt gene_snp_list.txt gene_snp_list_knn.txt; do
                    if [ -f ${RESULT_DIR}/chr${CHRM}/chunk000_${EXT} ]; then failedChunkReduce=true; echo "${TSHORT}: Error in chr${CHRM} ${EXT}"; fi
                done
            done
            #if [ "${failedChunkReduce}" = "false" ]; then echo "${TSHORT}: Chunk Reduce OK"; fi

            #failedPostResult=false
            #for FNAME in trans_eqtls.txt target_genes.txt; do
            #    if [ ! -f ${COLRES_DIR}/${FNAME} ]; then failedPostResult=true; fi
            #done
            #if [ "${failedPostResult}" = "false" ]; then echo "${TSHORT}: Post Result OK"; fi
            #if [ "${failedPostResult}" = "true" ];  then echo "${TSHORT}: Warning! No Result Created"; fi

        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
