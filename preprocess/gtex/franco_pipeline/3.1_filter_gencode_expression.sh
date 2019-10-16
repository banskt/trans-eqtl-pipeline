#!/bin/bash

source "./CONFIG.sh"

ALTTISSUEFILE="$HOME/Tejaas_Pipeline/devtools/tissues.table.txt"

while IFS='' read -r line || [ -n "$line" ]; do
	if [[ $line =~ ^[^\#] ]]; then
	    fullname=$(echo "$line" | cut -f 1 )
	    shortname=$(echo "$line" | cut -f 2 )
	    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
	    if [ $RUNPEER = true ]; then
	        for NCOV in `seq 0 5 $MAXNCOV`;
	        do
	            # get PEER covariates correcting for covariates above
	            PEER_EXPR_FILE="${COVDIR}/${shortname}_gtex.${NCOV}ncov_PEER_residuals.txt"
	            python filter_gencode_expr.py --gx ${PEER_EXPR_FILE} --donors ${DONORFILE} --dataset gtex
	        done
	    fi
	    EXPRFILE="${EXPROUTDIR}/gtex.normalized.expression.${shortname}.txt"
	    python filter_gencode_expr.py --gx ${EXPRFILE} --donors ${DONORFILE} --dataset gtex

	    LMEXPRFILE="${BASEDIR}/expression/norm_lmcorrected/gtex.normalized.expression.lmcorrected.${shortname}.txt"
	    python filter_gencode_expr.py --gx ${LMEXPRFILE} --donors ${DONORFILE} --dataset gtex

	    PEER1COVEXPR="${COVDIR}/${shortname}_gtex.1ncov_PEER_residuals.txt"
	    python filter_gencode_expr.py --gx ${PEER1COVEXPR} --donors ${DONORFILE} --dataset gtex
	fi

done < $TISSUEFILE
# done < $ALTTISSUEFILE

