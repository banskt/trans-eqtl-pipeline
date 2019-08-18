#!/bin/bash

#### variables used ####
# PYENV
# SRCTPM
# SRCREAD
# SRCPHENO 
# GENCODEFILE
#
# TFULL
# TSHORT
# TISSUEOUTDIR
# PREGXOUTDIR
# PREPROC_UTILSDIR
# COVARDIR
# NORMQCPY
#
# bCollectCovs
# bSelectTissue
# bNormalizeQC
#
# TPM_THRESHOLD
# COUNTS_THRESHOLD
# SAMPLE_FRAC_THRESHOLD
# QCMETHODS
# GXSELECTION 

## Collect all covariates. This is a placeholder. Change later
if [ "${bCollectCovs}" = "true" ]; then cp ${COVARDIR}/cov_gender_age_trischd.txt ${TISSUEOUTDIR}/covariates.txt; fi

## Select samples belonging to the particular tissue // Bash
if [ "${bSelectTissue}" = "true" ]; then 
    source ${PREPROC_UTILSDIR}/extract_selected_samples_tpm_counts
    ## column 15 (SMTSD) = tissue name, and column 29 (SMAFRZE) = RNASEQ
    tail -n +11 ${SRCPHENO} | awk -F $'\t' -v TISSUE="$TFULL" 'BEGIN {OFS = FS} $15 == TISSUE {print}' \
                            | awk -F $'\t' 'BEGIN {OFS = FS} $29 == "RNASEQ" {print $2}' \
                            > ${TISSUEOUTDIR}/selected_samples.txt
    extract_selected_samples_tpm_counts ${SRCTPM} ${SRCREAD} \
                            ${TISSUEOUTDIR}/selected_samples.txt ${PREGXOUTDIR}/rnaseq_samples.txt \
                            ${TISSUEOUTDIR}/all_genes_tpm.gct ${TISSUEOUTDIR}/all_genes_counts.gct
fi

## perform the QC and normalization // Python
if [ "${bNormalizeQC}" = "true" ]; then
    ${PYENV} ${NORMQCPY} --tpm ${TISSUEOUTDIR}/all_genes_tpm.gct \
                         --counts ${TISSUEOUTDIR}/all_genes_counts.gct \
                         --vcf_sample_list ${PREGXOUTDIR}/vcf_samples.list \
                         --out ${TISSUEOUTDIR}/gtexv8_${TSHORT}.txt \
                         --cov ${TISSUEOUTDIR}/covariates.txt \
                         --tpm_threshold ${TPM_THRESHOLD} \
                         --count_threshold ${COUNTS_THRESHOLD} \
                         --sample_frac_threshold ${SAMPLE_FRAC_THRESHOLD} \
                         --methods ${QCMETHODS} \
                         --gtf ${GENCODEFILE} \
                         --biotype ${GXSELECTION}
fi
