#!/bin/bash

source "./CONFIG.sh"

mkdir -p $LMOUTDIR

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    COVARS="${GTEXCOVDIR}/${shortname}_nopeer_covariates.txt"
    EXPRFILE="${EXPROUTDIR}/gtex.normalized.expression.${shortname}.txt"
    LMOUTFILE="${LMOUTDIR}/gtex.normalized.expression.lmcorrected.${shortname}.txt"
    if [ -e $EXPRFILE ] && [ -e $COVARS ]; then
        Rscript correct_expr_lm.R $EXPRFILE $LMOUTFILE $COVARS
    else
        echo "$fullname covariates not found"
    fi
done < $TISSUEFILE
