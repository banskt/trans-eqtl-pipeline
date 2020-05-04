#!/bin/bash

BIOTYPE_FILE="gencode.v26.gene_id_biotype_list.txt"
INFILE="$1"

while read LINE; do
  TISSUE=$( echo $LINE | awk '{print $1}' )
  VARIANT=$( echo $LINE | awk '{print $2}' )
  GENE=$(echo $LINE | awk '{print $4}')
  BIOTYPE=$( grep "${GENE}\." gencode.v26.gene_id_biotype_list.txt | awk '{print $2}' )
  echo -e "${TISSUE}\t${VARIANT}\t${GENE}\t${BIOTYPE}"
done < ${INFILE}
