#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./03_select_top_snps.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_top_snps


for EXPR_CORR in ${EXPRESSIONS}; do
    for MDATA in ${DATASETS}; do
        for KNN in ${KNNS}; do
            source ${DATALOAD}
            OUTDIR_DATA="${OUTDIR}/${MDATA}"

            if [ ! -z "$EXPRESSIONFILE" ]; then

                echo "Processing trans-eqtls for $MDATA"

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

                for NULL in ${TEJAAS_NULL}; do
                    if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                    if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi

                    for CUTOFF in ${TEJAAS_CUTOFF}; do
                        if [ ${DYNAMIC_SB} == "true" ]; then
                            for KEFF in ${KEFFS}; do
                                METHOD_VARIANT="${NULL}null_sbDynamic${KEFF}"
                                if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                                if [ "${bTejaas}" = "true" ];    then tejaas_top_snps "${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}" $CUTOFF; fi
                                if [ "${bTjsRandom}" = "true" ]; then tejaas_top_snps "${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}" $CUTOFF; fi
                            done
                        else
                            for SBETA in ${TEJAAS_SIGMA_BETA}; do
                                METHOD_VARIANT="${NULL}null_sb${SBETA}"
                                if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                                if [ "${bTejaas}" = "true" ];    then tejaas_top_snps "${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}" $CUTOFF; fi
                                if [ "${bTjsRandom}" = "true" ]; then tejaas_top_snps "${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}" $CUTOFF; fi
                            done
                        fi
                    done
                done
            fi
            unset_vars ${DATALOAD}
        done
    done
done

unset_vars ${EXTERNALLOAD}
unset_vars ${CONFIGFILE}
unset_vars PATHS
