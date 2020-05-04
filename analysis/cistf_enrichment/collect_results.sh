#!/bin/bash

CIS_INKB=$1

OUTDIR="${HOME}/trans-eQTL/dev-pipeline/analysis/distance_from_transcription_factor/${CIS_INKB}kb_window_tissue_specific"
RESDIR="${HOME}/trans_eqtl_results/gtex_v8_tejaas_permnull_sb0.1_knn"
COLFILE="${OUTDIR}/cistf_enrichment_${CIS_INKB}kb.txt"

echo -e "TISSUE\tN_TRANSEQTLS\tCISTF_FRAC\tENRICHMENT\tP_VALUE" > ${COLFILE}

for TISSUE in $( ls ${RESDIR}/ ); do

    RESFILE="${OUTDIR}/cistf_enrichment_${TISSUE}.txt"

    if [ -f ${RESFILE} ]; then
        NTRANS=$( wc -l ${RESDIR}/${TISSUE}/trans_eqtls.txt | cut -d" " -f1 )
        NTRANS=$(( NTRANS - 1 ))
        TFFRAC=$( cat ${RESFILE} | awk '{print $1}' )
        ENRICH=$( cat ${RESFILE} | awk '{print $2}' )
        PVALUE=$( cat ${RESFILE} | awk '{print $3}' )
        echo -e "${TISSUE}\t${NTRANS}\t${TFFRAC}\t${ENRICH}\t${PVALUE}" >> ${COLFILE}
    fi

done
