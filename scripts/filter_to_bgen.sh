#!/bin/bash

SRCDIR="/scratch/sbanerj/data/Cardiogenics"

for CHRM in {1..22}; do

    SRCFILE="${SRCDIR}/genotype/CG_${CHRM}.imputed.gz"
    SAMPLEFILE="${SRCDIR}/genotype/CG_plink.sample"
    IN_FILE="${SRCDIR}/genotype_qc/CG_dosages_filtered_${CHRM}.imputed.gz"
    SNPLIST="${SRCDIR}/genotype_qc/CG_filtered_snps_${CHRM}.txt"
    OUTFILE="${SRCDIR}/genotype_qc/CG_filtered_imputed_${CHRM}"

    zcat ${IN_FILE} | cut -d" " -f2 > ${SNPLIST}
    plink2 --gen ${SRCFILE} --sample ${SAMPLEFILE} --oxford-single-chr ${CHRM} --extract ${SNPLIST} --export bgen-1.2 --out ${OUTFILE}
done
