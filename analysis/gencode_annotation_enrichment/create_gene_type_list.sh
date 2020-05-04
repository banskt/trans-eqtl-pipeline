#!/bin/bash

GTFFILE=$1
OUTFILE=$2

if [ -n "${OUTFILE}" ]; then
  zcat $GTFFILE | awk 'BEGIN { FS=OFS="\t" } $3 == "gene" {print $9}' | \
                  awk 'BEGIN { FS=OFS="\t" } {split($1, meta, ";"); 
                                              split(meta[1], arr01, " "); 
                                              split(meta[2], arr02, " "); 
                                              gsub(/\"/, "", arr01[2]); 
                                              gsub(/\"/, "", arr02[2]); 
                                              gene_id = arr01[2];
                                              gene_type = arr02[2];
                                             } {print gene_id, gene_type}' \
                  > ${OUTFILE}
fi
