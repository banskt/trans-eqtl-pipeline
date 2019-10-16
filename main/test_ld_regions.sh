#!/bin/bash

LD_REGIONS_FILE="/cbscratch/franco/datasets/1KG_genomes/LD_regions2calculate.txt"

while IFS="" read -r line || [ -n "$p" ]
do
  # printf '%s\n' "$line"
  CHROM=`echo ${line} | cut -d" " -f 1`
  START=`echo ${line} | cut -d" " -f 2`
  END=`echo ${line} | cut -d" " -f 3`
  echo "CHROM ${CHROM} from ${START} to ${END}"
  
done < $LD_REGIONS_FILE