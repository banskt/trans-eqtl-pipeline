#!/bin/bash


$PYENV ${SCRIPTDIR}/filter_gencode_expr.py --gx ${EXPRFILE} --donors ${DONORFILE} --dataset gtex

if [ "${CORRECT_AGE}" = "true" ]; then
    $PYENV ${SCRIPTDIR}/filter_gencode_expr.py --gx ${LMOUTFILE_AGE} --donors ${DONORFILE} --dataset gtex
    PEERDIR="${PEEROUTDIR}_w_age"
    for NCOV in $NCOVS;
    do
        PEER_EXPR_FILE="${PEERDIR}/${TSHORT}_${NCOV}_PEER_residuals.txt"
        if [ -e $PEER_EXPR_FILE ]; then
            $PYENV ${SCRIPTDIR}/filter_gencode_expr.py --gx ${PEER_EXPR_FILE} --donors ${DONORFILE} --dataset gtex
        else
            "${PEER_EXPR_FILE} NOT FOUND"
        fi
    done
else
    $PYENV ${SCRIPTDIR}/filter_gencode_expr.py --gx ${LMOUTFILE} --donors ${DONORFILE} --dataset gtex
    PEERDIR="${PEEROUTDIR}"
    for NCOV in $NCOVS;
    do
        PEER_EXPR_FILE="${PEERDIR}/${TSHORT}_${NCOV}_PEER_residuals.txt"
        if [ -e $PEER_EXPR_FILE ]; then
            $PYENV ${SCRIPTDIR}/filter_gencode_expr.py --gx ${PEER_EXPR_FILE} --donors ${DONORFILE} --dataset gtex
        else
            "${PEER_EXPR_FILE} NOT FOUND"
        fi
    done
fi
