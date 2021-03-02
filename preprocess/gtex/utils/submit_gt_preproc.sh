#!/bin/bash

THISJOBDEPS="None"
SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/preprocess/gtex_${GTVERSION}/genotype/all_samples"
if [ ! -d ${SPECIFIC_JOBSUBDIR} ]; then mkdir -p ${SPECIFIC_JOBSUBDIR}; fi

for CHRM in {1..22}; do
    
    JOBNAME="vcf_split_filter_headerchange_annotation_chr${CHRM}_${RANDSTRING}"
    OUTFILE="${OUTFILE_BASE}_chr${CHRM}"
    ANNOTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot"
    
    if [ "${USE_SLURM}" = "true" ]; then
        __MASTERFILE="slurmfiles/vcf_split_filter_headerchange_annotation.sbatch"
        __SUBMITFILE="${SPECIFIC_JOBSUBDIR}/${JOBNAME}.sbatch"
    elif [ "${USE_LSF}" = "true" ]; then
        __MASTERFILE="lsffiles/vcf_split_filter_headerchange_annotation.sbatch"
        __SUBMITFILE="${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub"
    fi

    sed -e "s|_JOB_NAME|${JOBNAME}|g;
            s|_PYT_ENV_|${PYENV}|g;
            s|_VCFTOOLS|${VCFTOOLS}|g;
            s|_BGZIP___|${BGZIP}|g;
            s|_TABIX___|${TABIX}|g;
            s|_SRC_VCF_|${SRCVCF}|g;
            s|_ANNTFILE|${ANNOTFILE}|g;
            s|_OUT_FILE|${OUTFILE}|g;
            s|_ANNT_PY_|${CONVERT_ANNOT_PY}|g;
            s|_IMPT_PY_|${IMPUTE_MISSING_PY}|g;
            s|_VCFT_PY_|${VCF_FILTER_PY}|g;
            s|_CHRM_NUM|${CHRM}|g;
            s|_MAF_MIN_|${MAFMIN}|g;
            s|_RM_INDEL|${REMOVE_INDELS}|g;
            s|_RM_AMBIG|${REMOVE_AMBIGUOUS}|g;
            s|_IMP_MISS|${IMPUTE_MISSING}|g;
            s|_GTEX_VER_|${GTVERSION}|g;
           " ${__MASTERFILE} > ${__SUBMITFILE}

    if [ "${USE_SLURM}" = "true" ]; then
        GT_CHRM_JOBID=$( submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS} )
        JOBDEPS=$( add_deps "${JOBDEPS}" ${GT_CHRM_JOBID} )
    elif [ "${USE_LSF}" = "true" ]; then
        GT_CHRM_JOBID=$( submit_job_lsf ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS} )
        JOBDEPS=$( add_deps "${JOBDEPS}" ${GT_CHRM_JOBID} )
    fi

done
