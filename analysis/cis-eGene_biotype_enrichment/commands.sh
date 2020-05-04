#!/bin/bash

INFILE=$1
OUTPREFIX=$2

F_GTEX_UNIQ_BIOTYPES="gtex_v8/GTEx_Analysis_v8_uniq_cis-Genes.txt"
F_HUMANTF="../external/human_TF_annotation_gencode_v26.txt"

F_ALLTISSUE_BIOTYPE="cis-eGenes_tmplist_${OUTPREFIX}.txt"
F_ALLTISSUE_UNIQ_BIOTYPE="cis-eGenes_biotypes_${OUTPREFIX}.txt"
# F_TFCOUNTS="cis-eGenes_uniq_tfcounts_${OUTPREFIX}.txt"

./map_biotype.sh ${INFILE} > ${F_ALLTISSUE_BIOTYPE}
./calc_counts_biotypes.sh ${F_ALLTISSUE_BIOTYPE} ${F_HUMANTF} > ${F_ALLTISSUE_UNIQ_BIOTYPE}
# ./calc_counts_transcription_factors.sh ${F_ALLTISSUE_BIOTYPE} ${F_GTEX_UNIQ_BIOTYPES} ${F_HUMANTF} > ${F_TFCOUNTS}
