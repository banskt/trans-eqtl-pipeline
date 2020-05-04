#!/bin/bash

BEDFILE=$1
LABELFILE=$2

NEWFILE="$( basename $BEDFILE .bed )_relabel.bed"
RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

COMMAND=$( 
echo -ne "sed -e \""
while read LINE; do 
    label=$( echo ${LINE} | cut -d" " -f1 )
    index=$( echo ${LINE} | cut -d" " -f2 )
    echo -ne "s|${label}|${index}|g;"
done < ${LABELFILE}
echo "\"" )

echo ${COMMAND} ${BEDFILE} > sed_command_${RANDSTRING}.sh
chmod +x sed_command_${RANDSTRING}.sh
./sed_command_${RANDSTRING}.sh > ${NEWFILE}

rm -f sed_command_${RANDSTRING}.sh
