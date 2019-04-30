#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_process_chunks.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source EXTERNAL
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/gx_preproc_string
TEJAAS_FIXPVAL_PY=$(abs_path "${SCRIPTDIR}/tejaas_fixpvals.py")

function tejaas_fixpval() {
    local PYENV=$1
    local PYSCRIPT=$2
    local DIR=$3
    local CWD=$( pwd )
    cd $DIR
    #echo $PYENV ${PYSCRIPT} $DIR
    if [ -f rr.txt ]; then
        mv -v rr.txt rr_old.txt
    fi
    ${PYENV} ${PYSCRIPT} --infile $DIR/rr_old.txt --outfile $DIR/rr.txt
    cd $CWD
}

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
                echo "fixing chr${CHRM}"
                for NULL in ${TEJAAS_NULL}; do
                    if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                    if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                    for SBETA in ${TEJAAS_SIGMA_BETA}; do
                        if [ "${bTejaas}" = "true" ];    then tejaas_fixpval "${PYTHON36}" "${TEJAAS_FIXPVAL_PY}" "${OUTDIR_DATA}/tejaas/${NULL}null_sb${SBETA}/chr${CHRM}"; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_fixpval "${PYTHON36}" "${TEJAAS_FIXPVAL_PY}" "${OUTDIR_DATA}/tejaas_rand/${NULL}null_sb${SBETA}/chr${CHRM}"; fi
                    done
                done
            done
        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
