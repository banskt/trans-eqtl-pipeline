#!/bin/bash

GTEXDIR="$1"

RAND=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

for file in ${GTEXDIR}/GTEx_Analysis_v8_eQTL/*.egenes.txt.gz; do 
    zcat $file | cut -f1 | tail -n +2 >> tmp_${RAND}; 
    echo $file; 
done
sort tmp_${RAND} > tmp2_${RAND}
cat tmp2_${RAND} | uniq -c | awk '{printf "%s\t%d\n", $2, $1}' > tmp3_${RAND}
echo "Checking biotypes"
while read LINE; do
  GENE=$(echo $LINE | awk '{print $1}')
  NOCC=$(echo $LINE | awk '{print $2}')
  BIOTYPE=$( grep "${GENE}\>" ../gencode.v26.gene_id_biotype_list.txt | awk '{print $2}' )
  echo -e "${GENE}\t${NOCC}\t${BIOTYPE}"
done < tmp3_${RAND} > GTEx_Analysis_v8_uniq_cis-Genes.txt


