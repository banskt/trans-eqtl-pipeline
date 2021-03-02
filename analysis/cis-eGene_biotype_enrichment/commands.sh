#!/bin/bash

INFILE=$1
ALTSB_INFILE=$2
OUTPREFIX=$3
ALTSB_TISSUES="haa pan spl wb"

F_GTEX_UNIQ_BIOTYPES="gtex_v8/GTEx_Analysis_v8_uniq_cis-Genes.txt"
F_HUMANTF="../external/human_TF_annotation_gencode_v26.txt"
F_INFILE="ciseqtl_targets_pc_lncRNA_${OUTPREFIX}.txt"

F_ALLTISSUE_BIOTYPE="cis-eGenes_tmplist_${OUTPREFIX}.txt"
F_ALLTISSUE_UNIQ_BIOTYPE="cis-eGenes_biotypes_${OUTPREFIX}.txt"
# F_TFCOUNTS="cis-eGenes_uniq_tfcounts_${OUTPREFIX}.txt"

for TISSUE in ${ALTSB_TISSUES}; do
    if [ ! -z ${SEARCHSTR} ]; then SEARCHSTR="${SEARCHSTR}\|"; fi
    SEARCHSTR="${SEARCHSTR}${TISSUE}"
done
echo ${SEARCHSTR}
grep -vw "${SEARCHSTR}" ${INFILE} > ${F_INFILE}
grep -w "${SEARCHSTR}" ${ALTSB_INFILE} >> ${F_INFILE}

./map_biotype.sh ${F_INFILE} > ${F_ALLTISSUE_BIOTYPE}
./calc_counts_biotypes.sh ${F_ALLTISSUE_BIOTYPE} ${F_HUMANTF} > ${F_ALLTISSUE_UNIQ_BIOTYPE}
# rm -f ${F_INFILE}
# ./calc_counts_transcription_factors.sh ${F_ALLTISSUE_BIOTYPE} ${F_GTEX_UNIQ_BIOTYPES} ${F_HUMANTF} > ${F_TFCOUNTS}
