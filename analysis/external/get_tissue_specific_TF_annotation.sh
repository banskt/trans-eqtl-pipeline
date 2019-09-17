#!/bin/bash

TISSUEFILE="/usr/users/sbanerj/trans-eQTL/dev-pipeline/main/tissue_table.txt"
OUTDIR="tissue_specific_TFs"

cut -f1 human_TF_annotation_gencode_v26.txt | tail -n +2 > ${OUTDIR}/all_TF.txt

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        EXPRESSIONFILE="/scratch/sbanerj/trans-eqtl/input/gtex_v8/preprocess_gx/${TBASE}/gtexv8_${TSHORT}_raw_protein_coding.txt"
        echo ${TBASE}

        cut -f1 ${EXPRESSIONFILE} | tail -n +2 | cut -d"." -f1 > ${OUTDIR}/${TSHORT}_expressed.txt
        grep -Fwf ${OUTDIR}/all_TF.txt ${OUTDIR}/${TSHORT}_expressed.txt > ${OUTDIR}/${TSHORT}_TF_expressed.txt

        head -n 1 human_TF_annotation_gencode_v26.txt > ${OUTDIR}/human_TF_annotation_gencode_v26_${TSHORT}.txt
        grep -Fwf ${OUTDIR}/${TSHORT}_TF_expressed.txt human_TF_annotation_gencode_v26.txt >> ${OUTDIR}/human_TF_annotation_gencode_v26_${TSHORT}.txt

        rm -f ${OUTDIR}/${TSHORT}_expressed.txt
        rm -f ${OUTDIR}/${TSHORT}_TF_expressed.txt

    fi

done < ${TISSUEFILE}

rm -f ${OUTDIR}/all_TF.txt
