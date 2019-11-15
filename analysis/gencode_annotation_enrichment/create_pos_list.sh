#!/bin/bash

GTFFILE=$1
BIOTYPE=$2
OUTFILE=$3

if [ -n "${OUTFILE}" ]; then
    RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
    TMP_SEL_FILE="tmp_sel_${RANDSTRING}.txt"
    zcat $GTFFILE | awk 'BEGIN { FS=OFS="\t" } $3 == "gene"' | \
                    awk -v biotype="${BIOTYPE}" 'BEGIN { FS=OFS="\t" } 
                                                       {split($9, meta, ";");
                                                       for (i in meta) {
                                                         gsub (/^ */, "", meta[i]); 
                                                         if (meta[i] ~ /^gene_type/) {
                                                           split(meta[i], arr01, " "); 
                                                           gsub(/\"/, "", arr01[2]); 
                                                           gene_type = arr01[2];
                                                         }
                                                       }} 
                                                       gene_type == biotype {print $1, $4, $5, $9}' \
                     > ${TMP_SEL_FILE}
    echo -e 'ensembl_id\tchrom\tstart\tend\tname' > ${OUTFILE}
    awk -f reformat_gtf_selection.awk ${TMP_SEL_FILE} >> ${OUTFILE}
    rm -f ${TMP_SEL_FILE}
fi
## awk 'BEGIN { FS=OFS="\t" } $3 == "gene"' gencode.v26.annotation.gtf | cut -d";" -f2 | cut -d" " -f3 | sed "s/\"//g" | sort | uniq -c | sort -n
