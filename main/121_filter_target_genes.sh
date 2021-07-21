#!/bin/bash -e

module load intel/compiler/64/2020/19.1.2
module load intel/mkl/64/2020/2.254

echo "This script was never finished, because we don't really want to filter lead SNPs in the target gene list"
echo "If you want to finish it, call your script in line 65"
exit 0;

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
INFILES=""
for CUTOFF in ${TEJAAS_CUTOFF}; do
    for EXPR_CORR in ${EXPRESSIONS}; do
        for MDATA in ${DATASETS}; do
            # source ${DATALOAD}

            DATATYPE=`echo ${MDATA} | cut -d'-' -f1`
            if [ "$DATATYPE" == "fhs" ]; then
                TISSUEID="fhs"
            else
                TISSUEID=`echo ${MDATA} | cut -d'-' -f2`
            fi

            # This looks stupid, but is to be able to run LD across chromosomes always on the same dataset
            # If multiple datasets are to be run, then this will raise an error
            if [ "${CURRENT_DTYPE}" = "" ]; then CURRENT_DTYPE=$DATATYPE; fi
            if [ "${CURRENT_DTYPE}" != "${DATATYPE}" ]; then echo "Different datasets! abort"; exit; fi

            OUTDIR_DATA="${SOURCEDIR}/${EXPR_CORR}/summary_${CUTOFF}/${TISSUEID}"

            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for KNN in ${KNNS}; do

                    _VARIANT=""
                    if [ "${KNN}" -gt "0" ]; then _VARIANT="${_VARIANT}_knn${KNN}"; fi
                    if [ "${CISMASK}" != "true" ]; then _VARIANT="${_VARIANT}_nocismask"; fi
                    if [ "${CROSSMAP}" == "true" ]; then _VARIANT="${_VARIANT}_crossmap"; fi
                    if [ "${MAGIC_SQRT}" == "true" ]; then _VARIANT="${_VARIANT}_sqrt"; fi
                    if [ "${NOGTKNN}" == "true" ]; then _VARIANT="${_VARIANT}_nogtknn"; fi

                    for SBETA in $TEJAAS_SIGMA_BETA_PERM; do
                        METHOD_VARIANT="${NULL}null_sb${SBETA}"
                        if [ "${_VARIANT}" != "" ]; then METHOD_VARIANT="${METHOD_VARIANT}${_VARIANT}"; fi
                        if [ "${bTejaas}" = "true" ];   then
                            NEWFILE="${OUTDIR_DATA}/tejaas/${METHOD_VARIANT}/trans_eqtls.txt"
                            #### Run here the target gene filtering
                        fi
                        if [ "${bTjsRandom}" = "true" ];then
                            NEWFILE="${OUTDIR_DATA}/tejaas_rand/${METHOD_VARIANT}/trans_eqtls.txt" 
                            #### Run here the target gene filtering
                        fi
                    done
                done
            done
        done
        # LDFILE="${GENO_DIR}/ldmap/chr${CHRM}_${CURRENT_DTYPE}.geno.ld"
        # LDFILE="${GENO_DIR}/SHAPEIT2_ldmap_${LDWINDOW}_${LDMIN_R2}/chr{:d}_${CURRENT_DTYPE}.geno.ld"
        LDFILE="${GENO_DIR}/${LDFILE_DIR}/chr{:d}_${CURRENT_DTYPE}.geno.ld"
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
