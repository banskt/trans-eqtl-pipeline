#!/bin/bash

source "./CONFIG.sh"


# TISSUEFILE="tissues.table.txt"
# ENV=/home/fsimone/myenv/bin

# READSFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
# DONORFILE="../donor_ids.fam"
# EXPROUTDIR="normalized_expr"

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    echo "Processing Tissue: $fullname"
    RPKMINFILE="${RPKMOUTDIR}/GTEx_Data_20181031_RNAseq_RNASeQCv1.1.8_gene_rpkm_${shortname}.gct"
    $ENV/python do_expression_normalization.py --rpkm $RPKMINFILE --counts $READSFILE --tissue $shortname --donors $DONORFILE --outdir $EXPROUTDIR
done < $TISSUEFILE
