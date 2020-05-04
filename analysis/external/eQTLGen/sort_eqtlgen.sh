#!/bin/bash

INFILE=$1
OUTFILE=$2

head -n 1 ${INFILE} > ${OUTFILE}

for CHRM in $( seq 1 22 ); do
    echo $CHRM
    awk -v OFS='\t' -v CHRM=${CHRM} '$3 == CHRM {print $2, $3, $4, $5, $6, $8, $1}' ${INFILE} | sort -k3 -n >> ${OUTFILE}
done
