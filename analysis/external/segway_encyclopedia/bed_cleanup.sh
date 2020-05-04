#!/bin/bash

FILE=$1
NEWFILE="$( basename $FILE .bed )_clean.bed"

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
TMPDIR=bed_split_${RANDSTRING}
mkdir $TMPDIR
echo ${FILE}

for CHRM in {1..22}; do
    echo $CHRM
    awk -v CHRM=chr"$CHRM" '$1 == CHRM {print}' ${FILE} > ${TMPDIR}/chr${CHRM}.bed
done

if [ -f ${NEWFILE} ]; then rm -f ${NEWFILE}; fi
head -n 1 ${FILE} > ${NEWFILE}
for CHRM in {1..22}; do
    cat ${TMPDIR}/chr${CHRM}.bed >> ${NEWFILE}
done

rm -rf ${TMPDIR}
