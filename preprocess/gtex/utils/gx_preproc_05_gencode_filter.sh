#!/bin/bash

ALLGXFILES=()

EXPR_TYPES="qn tmm tpms rpkms"
CORR_TYPES="lasso cclm"

# normalized uncorrected
ALLGXFILES+=( )

for EXPR_TYPE in ${EXPR_TYPES}; do
    if [ "${EXPR_TYPE}" == "tpms" ] || [ "${EXPR_TYPE}" == "rpkms" ]; then
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${TSHORT}_${EXPR_TYPE}_qcfilter.txt" )
    else
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${TSHORT}_${EXPR_TYPE}.txt" )
    fi

    for CORR_TYPE in ${CORR_TYPES}; do        
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${TSHORT}_${EXPR_TYPE}_${CORR_TYPE}.txt" )
    done
done

# ALLGXFILES+=( "$LMOUTFILE_AGE" )
# for NPEER in ${GXNPEERS}; do
#     THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
#     PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
#     ALLGXFILES+=( "$PEER_EXPR_FILE" )
# done


## PEER correction on normalized gene expression (using other covariates)
# for bAgeCorrect in ${GXAGECORR}; do
#     if [ "${bAgeCorrect}" = "true" ]; then
#         ALLGXFILES+=( "$LMOUTFILE_AGE" )
#         for NPEER in ${GXNPEERS}; do
#             THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
#             PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
#             ALLGXFILES+=( "$PEER_EXPR_FILE" )
#         done
#     else
#         ALLGXFILES+=( "$LMOUTFILE" )
#         for NPEER in ${GXNPEERS}; do
#             THIS_PEEROUTDIR="${PEEROUTDIR}/covar"
#             PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
#             ALLGXFILES+=( "$PEER_EXPR_FILE" )

#             THIS_PEEROUTDIR="${PEEROUTDIR}/normalized"
#             PEER_EXPR_FILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
#             ALLGXFILES+=( "$PEER_EXPR_FILE" )
#         done
#     fi      
# done

for THIS_GXFILE in ${ALLGXFILES[@]}; do
    echo $THIS_GXFILE
    if [ -e $THIS_GXFILE ]; then 
        ${PYENV} ${GENCODEFILTERPY} --gx ${THIS_GXFILE} \
                                    --donors ${DONORFILE} --dataset gtex \
                                    --gtf ${GENCODEFILE} \
                                    --biotype ${GXSELECTION} # --out ${TISSUEGXFILE} 
    fi
done
