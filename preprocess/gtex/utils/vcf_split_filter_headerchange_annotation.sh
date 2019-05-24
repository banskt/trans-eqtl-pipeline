#!/bin/bash

THISJOBDEPS="None"
SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/preprocess/gtex/genotype/all_samples"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

for CHRM in {1..22}; do
    
    JOBNAME="gtex_vcf_split_filter_headerchange_annotation_chr${CHRM}_${RANDSTRING}"
    OUTFILE="${OUTFILE_BASE}_chr${CHRM}"
    ANNOTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot"
    
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
           " ${MASTER_BSUBDIR}/gtex_vcf_split_filter_headerchange_annotation.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

    echo "Job file: ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub"
    submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
    JOBDEPS=$( add_deps "${JOBDEPS}" ${JOBNAME} )
done
