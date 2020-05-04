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

            PP_JOBNAME="${TSHORT}_postprocess_${RANDSTRING}"
            RESULT_DIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}"
            COLRES_DIR="${TEQTL_OUTDIR}/${TSHORT}"
            CHRM_NTOT_FILE="$( dirname ${GENO_FMT} )/ntot_per_chromosome.txt"
            CHUNK_REDUCE_UTIL="${UTILSDIR}/tejaas_chunk_reduce"

            if [ ! -d ${SPECIFIC_JOBSUBDIR} ]; then mkdir -p ${SPECIFIC_JOBSUBDIR}; fi
            if [ ! -d ${COLRES_DIR} ]; then mkdir -p ${COLRES_DIR}; fi

            sed -e "s|_JOB_NAME|${PP_JOBNAME}|g;
                    s|_RES_DIR_|${RESULT_DIR}|g;
                    s|_COL_DIR_|${COLRES_DIR}|g;
                    s|_CHR_NTT_|${CHRM_NTOT_FILE}|g;
                    s|_TJC_RED_|${CHUNK_REDUCE_UTIL}|g;
                    s|_CHR_NUM_|\"${CHRNUMS}\"|g;
                    s|_TJ_PCUT_|${PP_TEJAAS_PVALCUT}|g;
                    s|_NMAX_JB_|${MAX_NSNP_PERJOB}|g;
                    s|_bCHNK_R_|${bChunkReduce}|g;
                    s|_bPST_RS_|${bPostResult}|g;
                   " ${MASTER_BSUBDIR}/postprocess.sbatch > ${SPECIFIC_JOBSUBDIR}/${PP_JOBNAME}.sbatch

            if [ ${PP_USE_SLURM} == "true" ]; then
                PP_JOBID=$( submit_job ${SPECIFIC_JOBSUBDIR} ${PP_JOBNAME} None )
                JOBDEPS=$( add_deps "${JOBDEPS}" ${PP_JOBID} )
                echo "${PP_JOBID}: Postprocess job for ${TBASE} (${PP_JOBNAME})"
            else
                echo "Postprocess for ${TBASE}."
                grep -v "#SBATCH" ${SPECIFIC_JOBSUBDIR}/${PP_JOBNAME}.sbatch > ${SPECIFIC_JOBSUBDIR}/${PP_JOBNAME}.sh
                source ${SPECIFIC_JOBSUBDIR}/${PP_JOBNAME}.sh
                rm -f ${SPECIFIC_JOBSUBDIR}/${PP_JOBNAME}.sbatch
            fi
        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
