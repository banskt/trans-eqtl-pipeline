#!/bin/bash

INFILE=$1
OUTFILE=$2

head -n 1 ${INFILE} > ${OUTFILE}

for CHRM in $( seq 1 22 ); do
    echo $CHRM
    awk -v OFS='\t' -v CHRM=${CHRM} '$2 == CHRM {print $1, $2, $3}' ${INFILE} | sort -k3 -n | uniq >> ${OUTFILE}
done
