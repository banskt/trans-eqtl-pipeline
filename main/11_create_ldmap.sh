#!/bin/bash

DATATYPE=$1

echo "Input data type: ${DATATYPE}"

if [ -z ${DATATYPE} ]; then
    echo "Fatal! No datatype specified";
    echo "Use this script as: ./11_create_ldmap.sh DATATYPE (Cardiogenics / GTEx)"
    exit 1
fi

if [ ! "${DATATYPE}" == "Cardiogenics" ]; then
    echo "Fatal! No such datatype.";
    echo "Please check your argument."
    echo "Use this script as: ./11_create_ldmap.sh DATATYPE (Cardiogenics / GTEx)"
    exit 1
fi

source PATHS
if [ "${DATATYPE}" == "Cardiogenics" ]; then
    BGEN_FMT=${CARDIO_BGEN_FMT}
elif [ "${DATATYPE}" == "Cardiogenics" ]; then
    BGEN_FMT=${GTEX_BGEN_FMT}
fi

OUTDIR_DATA="${OUTDIR}/ldmap/${DATATYPE}"
JOBDIR_DATA="${JOBSUBDIR}/ldmap/${DATATYPE}"
if [ ! -d ${OUTDIR_DATA} ]; then mkdir -p ${OUTDIR_DATA}; fi
if [ -d ${JOBDIR_DATA} ]; then rm -rf ${JOBDIR_DATA}; fi; mkdir -p ${JOBDIR_DATA}

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

cd ${JOBDIR_DATA}
for CHRM in 1 2 3 6; do

    JOBNAME="ld_${DATATYPE}_${CHRM}_${RANDSTRING}"
    BGENFILE="${BGEN_FMT/\[CHRM\]/${CHRM}}"
    BASEPREF=`basename ${BGENFILE} .bgen`
    CORRFILE="${OUTDIR_DATA}/${BASEPREF}.bcor"

    sed -e "s|_JOB_NAME|${JOBNAME}|g;
            s|_LD_STOR_|${LDSTORE}|g;
            s|_GT_FILE_|${BGENFILE}|g;
            s|_COR_FIL_|${CORRFILE}|g;
           " ${MASTER_BSUBDIR}/create_ldmap.bsub > ${JOBNAME}.bsub
    bsub < ${JOBNAME}.bsub

done
cd ${CURDIR}

source unset_variables.sh PATHS
