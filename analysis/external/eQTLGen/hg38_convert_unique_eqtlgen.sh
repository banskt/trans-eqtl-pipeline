#!/bin/bash

INFILE=$1
HG38FILE=$2
OUTFILE=$3

head -n 1 ${INFILE} > ${OUTFILE}
tail -n +2 ${INFILE} > tmp_hg19_noheader.txt

LINENUM=1
while read TEQTL; do
    NEWPOS=$( sed "${LINENUM}q;d" $HG38FILE | cut -f2 )
    echo ${TEQTL} | awk -v OFS='\t' -v MNEWPOS=${NEWPOS} '{print $1, $2, MNEWPOS}'
    LINENUM=$(( LINENUM + 1 ))
done < tmp_hg19_noheader.txt > tmp_hg38_noheader.txt

for CHRM in $( seq 1 22 ); do
    echo $CHRM
    awk -v OFS='\t' -v CHRM=${CHRM} '$2 == CHRM {print $1, $2, $3}' tmp_hg38_noheader.txt | sort -k3 -n | uniq >> ${OUTFILE}
done
rm -f tmp_hg19_noheader.txt
rm -f tmp_hg38_noheader.txt
