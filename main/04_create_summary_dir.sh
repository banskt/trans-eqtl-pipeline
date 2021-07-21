#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./04_create_summary_dir.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_top_snps

#dummy variable
if [ "$TISSUEFILE" == "fhs" ]; then
    MDATA="fhs"
elif [ "$TISSUEFILE" == "geu" ]; then
    MDATA="geu"
else
    MDATA="gtex_v8-no_file"
fi

for EXPR_CORR in ${EXPRESSIONS}; do
    PREPROCS=" "
    source ${DATALOAD}
    OUTDIR_DATA="${OUTDIR}"

    for KNN in ${KNNS}; do

        echo "Processing trans-eqtls results"

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

        if [ "${NOGTKNN}" == "true" ]; then
            _VARIANT="${_VARIANT}_nogtknn"
        fi

        for NULL in ${TEJAAS_NULL}; do
            if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
            if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi

            for CUTOFF in ${TEJAAS_CUTOFF}; do
                if [ ${DYNAMIC_SB} == "true" ]; then
                    for KEFF in ${KEFFS}; do
                        METHOD_VARIANT="${NULL}null_sbDynamic${KEFF}"
                        if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                        PREPROCS="${PREPROCS}${METHOD_VARIANT} "
                    done
                else
                    for SBETA in ${TEJAAS_SIGMA_BETA}; do
                        METHOD_VARIANT="${NULL}null_sb${SBETA}"
                        if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                        PREPROCS="${PREPROCS}${METHOD_VARIANT} "
                    done
                fi
            done
        done
    done

    echo $PREPROCS
    echo $TISSUEFILE
    echo $EXPR_CORR
    echo $OUTDIR_DATA
    tejaas_summary $OUTDIR_DATA "${PREPROCS}" $TISSUEFILE $DATATYPE 

    unset_vars ${DATALOAD}
done

unset_vars ${EXTERNALLOAD}
unset_vars ${CONFIGFILE}
unset_vars PATHS
