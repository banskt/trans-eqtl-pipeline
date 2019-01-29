#!/bin/bash

THISJOBDEPS="None"
SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/preprocess/gtex/all_samples"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

for CHRM in {1..22}; do
    
    JOBNAME="gtex_vcf_split_filter_headerchange_annotation_chr${CHRM}_${RANDSTRING}"
    OUTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_dbSNP135_maf${MAFMIN#*.}_chr${CHRM}"
    ANNOTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot"
    
    sed -e "s|_JOB_NAME|${JOBNAME}|g;
            s|_VCFTOOLS|${VCFTOOLS}|g;
            s|_BGZIP___|${BGZIP}|g;
            s|_TABIX___|${TABIX}|g;
            s|_SRC_VCF_|${SRCVCF}|g;
            s|_CHRM_NUM|${CHRM}|g;
            s|_OUT_FILE|${OUTFILE}|g;
            s|_MAF_MIN_|${MAFMIN}|g;
            s|_ANNOT_PY|${CONVERT_ANNOT}|g;
            s|_ANNTFILE|${ANNOTFILE}|g;
           " ${MASTER_BSUBDIR}/gtex_vcf_split_filter_headerchange_annotation.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

    echo "Job file: ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub"
    submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
    JOBDEPS=`add_deps "${JOBDEPS}" ${JOBNAME}`
done
