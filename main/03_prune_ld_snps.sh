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
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/unset_vars
source EXTERNAL

source ${DATALOAD}

INPUTFILES=""
LDOUTFILES=""
while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        COLRES_DIR="${TEQTL_OUTDIR}/${TSHORT}"
        TEQTLFILE="${COLRES_DIR}/trans_eqtls.txt"
        if [ -f ${TEQTLFILE} ]; then
            INPUTFILES="${INPUTFILES}${TEQTLFILE} "
            LDOUTFILES="${LDOUTFILES}${COLRES_DIR}/trans_eqtls_ldpruned.txt "
        fi

    fi
done < ${TISSUEFILE}

SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/postprocess"
LD_JOBNAME="ld_pruning_${RANDSTRING}"

sed -e "s|_JOB_NAME|${LD_JOBNAME}|g;
        s|_PYT_ENV_|${PYTHON36}|g;
        s|_LDX_FMT_|${LDFILE_FMT}|g;
        s|_LDPR_PY_|${LDPRUNE_PY}|g;
        s|_IN_FILS_|\"${INPUTFILES}\"|g;
        s|_OUT_FLS_|\"${LDOUTFILES}\"|g;
       " ${MASTER_BSUBDIR}/ldprune.sbatch > ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sbatch

if [ ${LD_USE_SLURM} == "true" ]; then
    LD_JOBID=$( submit_job ${SPECIFIC_JOBSUBDIR} ${LD_JOBNAME} None )
    JOBDEPS=$( add_deps "${JOBDEPS}" ${LD_JOBID} )
    echo "${LD_JOBID}: LD pruning for all tissues (${LD_JOBNAME})"
else
    grep -v "#SBATCH" ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sbatch > ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sh
    #source ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sh
    echo ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sh
    #rm -f ${SPECIFIC_JOBSUBDIR}/${LD_JOBNAME}.sbatch
fi

unset_vars ${CONFIGFILE}
unset_vars PATHS
