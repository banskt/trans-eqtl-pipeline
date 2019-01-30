#!/bin/bash
for CHRM in {1..22}; do
    OUTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}"
    vcftools --gzvcf ${SRCVCF} --chr ${CHRM} --recode --recode-INFO-all --out $OUTFILE
    sed "18q;d" ${OUTFILE}.recode.vcf | awk -F$'\t' '{for (i = 1; i <= NF; i++) {idmod=gensub(/^(GTEX-[^-]*).*/,"\\1","g", $i); print idmod} }' | awk 'BEGIN { ORS = "\t" } { print }' > newheader.txt
    bgzip -c ${OUTFILE}.recode.vcf > ${OUTFILE}.recode.vcf.gz
    tabix -r newheader.txt ${OUTFILE}.recode.vcf.gz > ${OUTFILE}.vcf.gz
    tabix -f ${OUTFILE}.vcf.gz
    rm -rf ${OUTFILE}.recode.vcf* newheader.txt
done
