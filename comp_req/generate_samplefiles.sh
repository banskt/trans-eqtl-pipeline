#!/bin/bash

SAMPFILE=$1
FILENAME=$(basename -- $SAMPFILE)
EXT="${FILENAME##*.}"
BASE="${FILENAME%.*}"

echo $FILENAME
echo $EXT
echo $BASE

NSAMPMAX="581" #as samples with genotype
NSAMPLES="200 300 400 500 ${NSAMPMAX}"

for NSAMP in $NSAMPLES; do
    echo $NSAMP;
    head -n $(( $NSAMP + 2 )) "${SAMPFILE}" > "samples/${BASE}_N${NSAMP}.txt"; 
done;
