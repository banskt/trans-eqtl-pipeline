#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_process_chunks.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_chunk_reduce

# for saikat
# source ${UTILSDIR}/gx_preproc_string
# source ${UTILSDIR}/tejaas_chunk_reduce.new

for EXPR_CORR in ${EXPRESSIONS}; do
    for MDATA in ${DATASETS}; do
        for KNN in ${KNNS}; do
            # source DATA
            source ${DATALOAD}
            OUTDIR_DATA="${OUTDIR}/${MDATA}"

            echo $EXPRESSIONFILE

            if [ ! -z "$EXPRESSIONFILE" ]; then

                echo "Processing chunks for $MDATA"

                _VARIANT=""
                if [ "${KNN}" -gt "0" ]; then
                    _VARIANT="${_VARIANT}_knn${KNN}"
                fi

                if [ "${CISMASK}" != "true" ]; then
                    _VARIANT="${_VARIANT}_nocismask"
                fi

                if [ "${CROSSMAP}" == "true" ]; then
                    _VARIANT="${_VARIANT}_crossmap"
                fi

                if [ "${MAGIC_SQRT}" == "true" ]; then
                    _VARIANT="${_VARIANT}_sqrt"
                fi

                if [ "${SHUFFLE_SPECIAL}" == "true" ]; then
                    _VARIANT="${_VARIANT}_shuf"
                fi

                if [ "${NOGTKNN}" == "true" ]; then
                    _VARIANT="${_VARIANT}_nogtknn"
                fi

                for CHRM in ${CHRNUMS}; do
                    if [ "${bTejaasJPA}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/jpa/chr${CHRM}"; fi
                    if [ "${bJPARandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/jpa/chr${CHRM}"; fi
                    for NULL in ${TEJAAS_NULL}; do
                        if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                        if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi

                        if [ ${DYNAMIC_SB} == "true" ]; then
                            for KEFF in ${KEFFS}; do
                                METHOD_VARIANT="${NULL}null_sbDynamic${KEFF}"
                                if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                                if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}/chr${CHRM}"; fi
                                if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}/chr${CHRM}"; fi
                            done
                        else
                            for SBETA in ${TEJAAS_SIGMA_BETA}; do
                                METHOD_VARIANT="${NULL}null_sb${SBETA}"
                                if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                                if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}/chr${CHRM}"; fi
                                if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}/chr${CHRM}"; fi
                                # if [ "${bTjsRandomN}" = "true" ]; then 
                                #     for r in `seq 1 $NRANDOMS`; do
                                #         tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand_${r}/${NULL}null_sb${SBETA}/chr${CHRM}"
                                #     done;
                                # fi
                                # if [ "${bTejaasPartition}" = "true" ]; then
                                #     for r in `seq 1 $REPLICAS`; do
                                #         tejaas_chunk_reduce "${OUTDIR_DATA}/${r}/tejaas_part1/${NULL}null_sb${SBETA}/chr${CHRM}"
                                #         tejaas_chunk_reduce "${OUTDIR_DATA}/${r}/tejaas_part2/${NULL}null_sb${SBETA}/chr${CHRM}"
                                #     done;
                                # fi
                            done
                        fi
                    done
                done
            fi
        done
    done
done
unset_vars PATHS
unset_vars ${EXTERNALLOAD}
unset_vars ${CONFIGFILE}

