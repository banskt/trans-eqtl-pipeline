#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi


source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/unset_vars
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THISJOBDEPS="None"

CURRENT_DTYPE=""
for CHRM in ${CHRNUMS}; do
    INFILES=""
    for EXPR_CORR in ${EXPRESSIONS}; do
        for MDATA in ${DATASETS}; do
            source ${DATALOAD}

            # This looks stupid, but is to be able to run LD across chromosomes always on the same dataset
            # If multiple datasets are to be run, then this will raise an error
            if [ "${CURRENT_DTYPE}" = "" ]; then CURRENT_DTYPE=$DATATYPE; fi
            if [ "${CURRENT_DTYPE}" != "${DATATYPE}" ]; then echo "Different datasets! abort"; exit; fi

            JOBDIR_DATA="${JOBSUBDIR}/ldprune/${DATATYPE}"
            if [ -d ${JOBDIR_DATA} ]; then rm -rf ${JOBDIR_DATA}; fi; mkdir -p ${JOBDIR_DATA}
            JOBNAME="ldprune_${DATATYPE}_${CHRM}_${RANDSTRING}"

            OUTDIR_DATA="${OUTDIR}/${MDATA}"
            # if [ "${bMatrixEqtl}" = "true" ];  then source ${UTILSDIR}/matrix_eqtl; fi
            # if [ "${bMEqtlRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; fi
            for SBETA in $TEJAAS_SIGMA_BETA_PERM; do
                if [ "${bTejaas}" = "true" ];   then
                    NEWFILE="${OUTDIR_DATA}/tejaas/permnull_sb${SBETA}/chr${CHRM}/rr.txt"
                    INFILES="${INFILES}${NEWFILE} " ; fi
                if [ "${bTjsRandom}" = "true" ];then
                    NEWFILE="${OUTDIR_DATA}/tejaas_rand/permnull_sb${SBETA}/chr${CHRM}/rr.txt" 
                    INFILES="${INFILES}${NEWFILE} " ; fi
            done
            # if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi
            # if [ "${bJPARandom}" = "true" ];   then SHUFFLE=true; RUNJPA=true; source ${UTILSDIR}/tejaas; fi                
            
        done
        LDFILE="${GENO_DIR}/ldmap/chr${CHRM}_${CURRENT_DTYPE}.geno.ld"
        unset_vars ${DATALOAD}
    done
    PYPRUNELD="${SCRIPTDIR}/prune_ld_snps.py"

    sed -e "s@_JOB_NAME@${JOBNAME}@g;
            s@_PYENV_@${PYENV}@g;
            s@_PYPRUNE_@${PYPRUNELD}@g;
            s@_INFILES_@${INFILES}@g;
            s@_CHRM_@${CHRM}@g;
            s@_LDFILE_@${LDFILE}@g;
            " ${MASTER_BSUBDIR}/prune_ld.bsub > ${JOBDIR_DATA}/${JOBNAME}.bsub

    submit_job ${JOBDIR_DATA} ${JOBNAME} ${THISJOBDEPS}
done



