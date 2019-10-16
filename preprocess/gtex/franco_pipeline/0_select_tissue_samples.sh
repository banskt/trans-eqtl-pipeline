#!/bin/bash

source "./CONFIG.sh"

# TISSUEFILE="tissues.table.txt"
# ENV=/home/fsimone/myenv/bin
# RPKMFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
# PHENOFILE="phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt"
# RPKMOUTDIR="rpkms"

mkdir -p $RPKMOUTDIR;
rm -f 00_commands.sh

while IFS='' read -r line || [ -n "$line" ]; do
    # echo $line
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    # echo "Processing Tisuee: $fullname"
    RPKMOUTFILE="${RPKMOUTDIR}/GTEx_Data_20181031_RNAseq_RNASeQCv1.1.8_gene_rpkm_${shortname}.gct"
    echo $ENV/python select_samples_from_tissue.py --input $RPKMFILE --output $RPKMOUTFILE --tissue=\"$fullname\" --pheno $SAMPLE_PHENOFILE >> 00_commands.sh
done < $TISSUEFILE

sh 00_commands.sh