#!/bin/bash

source PATHS
TISSUE_FILE="../main/tissue_table.txt"
BASEDIR="/cbscratch/franco/trans-eqtl"
PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
PYPRUNELD="${SCRIPTDIR}/prune_ld_snps_gw.py"

CONFIGS="protein_coding_lncRNA_gamma01_knn30_cut5e-8_fst_high protein_coding_lncRNA_gamma0006_knn30_cut5e-8_fst_high"
# CONFIGS="protein_coding_lncRNA_gamma0006_knn30_cut5e-8_fst_high"
LDFILE="/cbscratch/franco/datasets/gtex_v8/genotypes/SHAPEIT2_ldmap_200000_0.5/chr{:d}_gtex_v8.geno.ld"


TSHORTS=""
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
        TSHORTS="${TSHORTS} ${TSHORT}"
    fi
done < ${TISSUE_FILE}


for CONFIG in $CONFIGS; do
    echo $CONFIG
    INFILES=""
    RESDIR="${BASEDIR}/${CONFIG}"
    for TISSUE in $TSHORTS; do
        NEWFILE="${RESDIR}/${TISSUE}/trans_eqtls.txt" 
        INFILES="${INFILES}${NEWFILE} "
    done

    $PYENV $PYPRUNELD --infiles $INFILES --ldfile $LDFILE 
done




