#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS

OPTIM_TISSUES="haa spl pan wb"
declare -A map
for name in $OPTIM_TISSUES; do
  map["$name"]=1
done

##################
DEV="false"
##################

for EXPR_CORR in ${EXPRESSIONS}; do
    source ${DATALOAD}
    OUTDIR_WEB="$OUTDIR/summary_web"
    if [ "$DEV" = "true" ]; then
        OUTDIR_WEB="$OUTDIR/summary_dev"
    fi
    mkdir -p $OUTDIR_WEB

    for TSHORT in ${TSHORTS}; do
        mkdir -p "${OUTDIR_WEB}/${TSHORT}"

        if [[ ${map["$TSHORT"]} ]] ; then
            TEJAAS_SIGMA_BETA_PERM="0.006"
        else
            TEJAAS_SIGMA_BETA_PERM="0.1"
        fi

        echo "$TSHORT permnull_sb${TEJAAS_SIGMA_BETA}"
        cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/trans_eqtls.txt" "${OUTDIR_WEB}/${TSHORT}/"
        cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/target_genes_FDR.txt" "${OUTDIR_WEB}/${TSHORT}/target_genes.txt"
        # cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/target_genes.txt" "${OUTDIR_WEB}/${TSHORT}/"
        cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/trans_eqtls_ldpruned.txt" "${OUTDIR_WEB}/${TSHORT}/"
        cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/ld_regions.txt" "${OUTDIR_WEB}/${TSHORT}/"
        if [ "$DEV" = "true" ]; then
            cp -v "$OUTDIR/summary_5e-08/$TSHORT/tejaas/permnull_sb${TEJAAS_SIGMA_BETA_PERM}_knn30/snps_list.txt" "${OUTDIR_WEB}/${TSHORT}/"
        fi
    done;
done;