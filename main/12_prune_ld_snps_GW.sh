#!/bin/bash


module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174

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
source ${UTILSDIR}/unset_vars
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THISJOBDEPS="None"

CURRENT_DTYPE=""
for CUTOFF in ${TEJAAS_CUTOFF}; do
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

            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for KNN in ${KNNS}; do

                    _VARIANT=""
                    if [ "${KNN}" -gt "0" ]; then _VARIANT="${_VARIANT}_knn${KNN}"; fi
                    if [ "${CISMASK}" != "true" ]; then _VARIANT="${_VARIANT}_nocismask"; fi
                    if [ "${MAGIC_SQRT}" == "true" ]; then _VARIANT="${_VARIANT}_sqrt"; fi

                    if [ ${DYNAMIC_SB} == "true" ]; then
                        METHOD_VARIANT="${NULL}null_sbDynamic"
                        if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                        if [ "${bTejaas}" = "true" ];   then
                            NEWFILE="${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}/trans_eqtls_${CUTOFF}.txt"
                            INFILES="${INFILES}${NEWFILE} " ; fi
                        if [ "${bTjsRandom}" = "true" ];then
                            NEWFILE="${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}/trans_eqtls_${CUTOFF}.txt" 
                            INFILES="${INFILES}${NEWFILE} " ; fi
                    else
                        for SBETA in $TEJAAS_SIGMA_BETA_PERM; do
                            METHOD_VARIANT="${NULL}null_sb${SBETA}"
                            if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                            if [ "${bTejaas}" = "true" ];   then
                                NEWFILE="${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}/trans_eqtls_${CUTOFF}.txt"
                                INFILES="${INFILES}${NEWFILE} " ; fi
                            if [ "${bTjsRandom}" = "true" ];then
                                NEWFILE="${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}/trans_eqtls_${CUTOFF}.txt" 
                                INFILES="${INFILES}${NEWFILE} " ; fi
                        done
                    fi
                    # if [ "${bTejaasJPA}" = "true" ];   then RUNJPA=true; source ${UTILSDIR}/tejaas; fi
                    # if [ "${bJPARandom}" = "true" ];   then SHUFFLE=true; RUNJPA=true; source ${UTILSDIR}/tejaas; fi                
                done
            done
        done
        # LDFILE="${GENO_DIR}/ldmap/chr${CHRM}_${CURRENT_DTYPE}.geno.ld"
        LDFILE="${GENO_DIR}/ldmap_${LDWINDOW}_${LDMIN_R2}/chr{:d}_${CURRENT_DTYPE}.geno.ld"
        unset_vars ${DATALOAD}
    done
    PYPRUNELD="${SCRIPTDIR}/prune_ld_snps_gw.py"

    $PYENV $PYPRUNELD --infiles $INFILES --ldfile $LDFILE 

    # sed -e "s@_JOB_NAME@${JOBNAME}@g;
    #         s@_PYENV_@${PYENV}@g;
    #         s@_PYPRUNE_@${PYPRUNELD}@g;
    #         s@_INFILES_@${INFILES}@g;
    #         s@_CHRM_@${CHRM}@g;
    #         s@_LDFILE_@${LDFILE}@g;
    #         " ${MASTER_BSUBDIR}/prune_ld.bsub > ${JOBDIR_DATA}/${JOBNAME}.bsub

    # submit_job ${JOBDIR_DATA} ${JOBNAME} ${THISJOBDEPS}
done



