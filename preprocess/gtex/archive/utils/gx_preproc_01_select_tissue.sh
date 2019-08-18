#!/bin/bash

#### variables used ####
# SRCTPM
# SRCREAD
# SRCPHENO 
# TFULL
# TISSUEOUTDIR
# PREGXOUTDIR

## column 15 (SMTSD) = tissue name, and column 29 (SMAFRZE) = RNASEQ

tail -n +11 ${SRCPHENO} | awk -F $'\t' -v TISSUE="$TFULL" 'BEGIN {OFS = FS} $15 == TISSUE {print}' | awk -F $'\t' 'BEGIN {OFS = FS} $29 == "RNASEQ" {print $2}' > ${TISSUEOUTDIR}/selected_samples.txt
grep -n -Fwf ${TISSUEOUTDIR}/selected_samples.txt ${PREGXOUTDIR}/rnaseq_samples.txt | cut -d':' -f1 > ${TISSUEOUTDIR}/tmp.columns
TRGTFIELDS="1,2,$(tr '\n' ',' < ${TISSUEOUTDIR}/tmp.columns)"
TRGTFIELDS=${TRGTFIELDS::-1}
zcat ${SRCTPM} | head -n 2  > ${TISSUEOUTDIR}/all_genes_tpm.gct
zcat ${SRCREAD} | head -n 2  > ${TISSUEOUTDIR}/all_genes_counts.gct
zcat ${SRCTPM} | tail -n +3 | cut -f$TRGTFIELDS >> ${TISSUEOUTDIR}/all_genes_tpm.gct
zcat ${SRCREAD} | tail -n +3 | cut -f$TRGTFIELDS >> ${TISSUEOUTDIR}/all_genes_counts.gct
rm -f ${TISSUEOUTDIR}/tmp.columns

