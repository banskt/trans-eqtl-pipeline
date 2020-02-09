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
# GXOUTDIR
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

## Select samples belonging to the particular tissue // Bash
if [ "${bSelectTissue}" = "true" ]; then 
    source ${PREPROC_UTILSDIR}/extract_selected_samples_tpm_counts

    if [[ ${SRCPHENO} == *"phs000424.v6"* ]]; then
        ## GTEx v6 | column 15 (SMTSD) = tissue name, and column 28 (SMAFRZE) != FLAGGED
        tail -n +11 ${SRCPHENO} | awk -F $'\t' -v TISSUE="$TFULL" 'BEGIN {OFS = FS} $15 == TISSUE {print}' \
                                | awk -F $'\t' 'BEGIN {OFS = FS} $28 != "FLAGGED" {print $2}' \
                                > ${TISSUEOUTDIR}/selected_samples.txt

    elif [[ ${SRCPHENO} == *"phs000424.v8"* ]]; then
        ## GTEx v8 | column 15 (SMTSD) = tissue name, and column 29 (SMAFRZE) = RNASEQ
        tail -n +11 ${SRCPHENO} | awk -F $'\t' -v TISSUE="$TFULL" 'BEGIN {OFS = FS} $15 == TISSUE {print}' \
                                | awk -F $'\t' 'BEGIN {OFS = FS} $29 == "RNASEQ" {print $2}' \
                                > ${TISSUEOUTDIR}/selected_samples.txt
    fi

    extract_selected_samples_tpm_counts ${SRCTPM} ${SRCREAD} \
                            ${TISSUEOUTDIR}/selected_samples.txt ${PREGXOUTDIR}/rnaseq_samples.txt \
                            ${TISSUEOUTDIR}/all_genes_tpm.gct ${TISSUEOUTDIR}/all_genes_counts.gct
fi

## perform the QC and normalization // Python
if [ "${bNormalizeQC}" = "true" ]; then
    ${PYENV} ${NORMQCPY} --tpm ${TISSUEOUTDIR}/all_genes_tpm.gct \
                         --counts ${TISSUEOUTDIR}/all_genes_counts.gct \
                         --vcf_sample_list ${PREGXOUTDIR}/vcf_samples.list \
                         --out ${TISSUEOUTDIR}/gtex_${TSHORT}.txt \
                         --cov ${TISSUEOUTDIR}/gtex_covariates.txt \
                         --tpm_threshold ${TPM_THRESHOLD} \
                         --count_threshold ${COUNTS_THRESHOLD} \
                         --sample_frac_threshold ${SAMPLE_FRAC_THRESHOLD} \
                         --methods ${QCMETHODS} \
                         --gtf ${GENCODEFILE} \
                         --biotype ${GXSELECTION}
fi

## copy files to the expression folder
cp ${TISSUEOUTDIR}/gtex_${TSHORT}* ${GXOUTDIR}/
