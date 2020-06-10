#!/bin/bash

EXPRFILE=$1
FILENAME=$(basename -- $EXPRFILE)
EXT="${FILENAME##*.}"
BASE="${FILENAME%.*}"

echo $FILENAME
echo $EXT
echo $BASE

NGENEMAX=15673 #(+1, for the header)
NGENES="5000 10000 ${NGENEMAX}"

for NGENE in $NGENES; do
    echo $NGENE;
    head -n $(( $NGENE + 1 )) "${EXPRFILE}" > "expr/${BASE}_N${NGENE}.txt"; 
done;
