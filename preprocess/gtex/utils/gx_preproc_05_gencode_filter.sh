#!/bin/bash

ALLGXFILES=()

# normalized uncorrected
ALLGXFILES+=( "$NORMOUTFILE" )

## PEER correction on normalized gene expression (using other covariates)
for bAgeCorrect in ${GXAGECORR}; do
    if [ "${bAgeCorrect}" = "true" ]; then
        ALLGXFILES+=( "$LMOUTFILE_AGE" )
        for NPEER in ${GXNPEERS}; do
            THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
            PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
            ALLGXFILES+=( "$PEER_EXPR_FILE" )
        done
    else
        ALLGXFILES+=( "$LMOUTFILE" )
        for NPEER in ${GXNPEERS}; do
            THIS_PEEROUTDIR="${PEEROUTDIR}/covar"
            PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
            ALLGXFILES+=( "$PEER_EXPR_FILE" )

            THIS_PEEROUTDIR="${PEEROUTDIR}/normalized"
            PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
            ALLGXFILES+=( "$PEER_EXPR_FILE" )
        done
    fi      
done

for THIS_GXFILE in ${ALLGXFILES[@]}; do
    if [ -e $THIS_GXFILE ]; then 
        ${PYENV} ${GENCODEFILTERPY} --gx ${THIS_GXFILE} \
                                    --donors ${DONORFILE} --dataset gtex \
                                    --gtf ${GENCODEFILE} \
                                    --biotype ${GXSELECTION} # --out ${TISSUEGXFILE} 
    fi
done
