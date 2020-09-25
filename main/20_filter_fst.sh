#!/bin/bash

source PATHS
TISSUE_FILE="../main/tissue_table.txt"
BASEDIR="/cbscratch/franco/trans-eqtl"
FSTFILE="/cbscratch/franco/from_saikat/EUR-AFR.weir.fst"
PYFILTERFST="${SCRIPTDIR}/filter_fst_snps_gw.py"
PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"

CONFIGS="protein_coding_lncRNA_gamma01_knn30_cut5e-8 protein_coding_lncRNA_gamma0006_knn30_cut5e-8"

for CONFIG in $CONFIGS; do
    echo $CONFIG
    OUTDIR="${BASEDIR}/${CONFIG}_fst_high"
    RESDIR="${BASEDIR}/${CONFIG}"
    $PYENV $PYFILTERFST  --tissues ${TISSUE_FILE} \
                        --fst ${FSTFILE}\
                        --basedir ${RESDIR}\
                        --outdir ${OUTDIR}
done

