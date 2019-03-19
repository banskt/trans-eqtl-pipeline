#!/bin/bash -e

mkdir -p $LMOUTDIR;

if [ -e $EXPRFILE ] && [ -e $COVARS ] && [ -e $COVARS_AGE ]; then
    if [ "${CORRECT_AGE}" = "true" ]; then
        Rscript ${SCRIPTDIR}/correct_expr_lm.R $EXPRFILE $LMOUTFILE_AGE $COVARS_AGE
    else
        Rscript ${SCRIPTDIR}/correct_expr_lm.R $EXPRFILE $LMOUTFILE $COVARS
    fi
else
    echo "$TFULL covariates or $EXPRFILE not found"
    exit
fi