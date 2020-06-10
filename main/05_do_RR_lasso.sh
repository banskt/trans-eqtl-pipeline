#!/bin/bash

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./05_do_RR_lasso.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/submit_job

INPUTDIR_SUMMARY="/cbscratch/franco/trans-eqtl/protein_coding_lncRNA_gamma01_knn30_cut5e-8"

# for t in `grep -v "#" ../../main/tissues.txt | cut -f 2`; 
# do sbatch -p hh --exclusive -N 1 -t 6-00:00:00 -o lasso_${t}.out -e lasso_${t}.err ç
# --wrap="$HOME/opt/miniconda/3/envs/env3.6/bin/python RR_lasso.py 
                # --tissue ${t} 
                # --outdir /cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze/lasso_targets"; done

PYRRLASSO="${SCRIPTDIR}/RR_lasso.py"

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THISJOBDEPS="None"


for EXPR_CORR in ${EXPRESSIONS}; do
    for TSHORT in ${TSHORTS}; do
        MDATA="gtex_v8-${TSHORT}"
        for KNN in ${KNNS}; do
            source ${DATALOAD}
            JOBSUBDIR_DATA="${JOBSUBDIR}/lasso_targets"
            OUTDIR_DATA="${INPUTDIR_SUMMARY}/lasso_targets"

            if [ ! -d ${OUTDIR_DATA} ]; then mkdir -p ${OUTDIR_DATA}; fi;
            if [ ! -d ${JOBSUBDIR_DATA} ]; then mkdir -p ${JOBSUBDIR_DATA}; fi

            if [ ! -z "$EXPRESSIONFILE" ]; then

                echo "Submitting Lasso RR for $TSHORT"
                
                for CHRM in ${CHRNUMS}; do
                    # GENOTYPEFILE=${GENO_FMT/\[CHRM\]/${CHRM}}
                    GENOTYPEFILE=${GENO_VCF_FMT/\[CHRM\]/${CHRM}}
                    if [ "${DATATYPE}" == "fhs" ]; then
                        GENOTYPEFILE=${GENO_FMT/\[CHRM\]/${CHRM}}
                    fi

                    SPECIFIC_OUTDIR="${OUTDIR_DATA}/${TSHORT}"
                    SPECIFIC_JOBSUBDIR="${JOBSUBDIR_DATA}/${TSHORT}/chr${CHRM}"
                    JOBNAME="rr_lasso_${TSHORT}_chr${CHRM}_${RANDSTRING}"

                    if [ ! -d ${SPECIFIC_OUTDIR} ]; then mkdir -p ${SPECIFIC_OUTDIR}; fi;
                    if [ ! -d ${SPECIFIC_JOBSUBDIR} ]; then mkdir -p ${SPECIFIC_JOBSUBDIR}; fi

                    sed "s|_JOB_NAME|${JOBNAME}|g;
                         s|_PY_ENV_|${PYENV}|g;
                         s|_RLASSO_|${PYRRLASSO}|g;
                         s|_SAMFIL_|${SAMPLEFILE}|g;
                         s|_GXFILE_|${EXPRESSIONFILE}|g;
                         s|_GXCORR_|${EXPRCORRFILE}|g;
                         s|_GENINF_|${GENEINFOFILE}|g;
                         s|_VCFFIL_|${GENOTYPEFILE}|g;
                         s|_OUTDIR_|${SPECIFIC_OUTDIR}|g;
                         s|_BASEDR_|${INPUTDIR_SUMMARY}|g;
                         s|_CHRM_N_|${CHRM}|g;
                         s|_TISSUE_|${TSHORT}|g;
                         s|_KNN_K_|${KNN}|g;
                         " ${MASTER_BSUBDIR}/tejaas_lasso_targets.slurm > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.slurm

                    submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}

                done
            fi
            unset_vars ${DATALOAD}
        done
    done
done

unset_vars ${EXTERNALLOAD}
unset_vars ${CONFIGFILE}
unset_vars PATHS
