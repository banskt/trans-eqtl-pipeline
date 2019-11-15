#!/bin/bash

NTARGET=$1
SNPPROP=$2
OUTDIR=$3

if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi
SNPCOUNTS="${OUTDIR}/snp_per_chrom.txt"
if [ -f ${SNPCOUNTS} ]; then rm -f ${SNPCOUNTS}; fi

## How many SNPs in each chromosome?
NRANDS=()
i=0
for CHRM in {1..22}; do
    NRANDS[i]=$( sed "${CHRM}q;d" ${SNPPROP} | awk '{print $2}' )
    i=$(( i + 1 ))
done
NRANDSTOT=0
for i in "${!NRANDS[@]}"; do
    NRANDSTOT=$(( NRANDSTOT + NRANDS[$i] ))
done
for i in "${!NRANDS[@]}"; do
    NSNP=$( echo ${NRANDS[$i]} ${NTARGET} | awk -v NTOT="${NRANDSTOT}" 'M = ($1 / NTOT) * $2 + 0.5 {printf "%d\n", M}' )
    echo $NSNP >> ${SNPCOUNTS}
done

## Create the list of SNPs    
for CHRM in {1..22}; do
    INFILE="/usr/users/sbanerj/gtex_v8/genotype/all_samples/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_withDS_chr${CHRM}.snplist"
    OUTFILE="${OUTDIR}/chr${CHRM}.txt"
    NSNP=$( sed "${CHRM}q;d" ${SNPCOUNTS} )
    cat ${INFILE} | shuf | head -n ${NSNP} | cut -d"_" -f2 > ${OUTFILE}
done
