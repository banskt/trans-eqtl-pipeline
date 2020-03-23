#!/bin/bash -e

ALLGXFILES=()

EXPR_TYPES="tmm tpms"
CORR_TYPES="cclm"

# EXPR_TYPES="tpms"
# CORR_TYPES="cclm cclm_nopc"


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
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${TSHORT}_${EXPR_TYPE}_${CORR_TYPE}_peer${GXNPEERS}_PEER_residuals.txt" )
    done
done

for THIS_GXFILE in ${ALLGXFILES[@]}; do
    echo $THIS_GXFILE
    if [ -e $THIS_GXFILE ]; then 
        ${PYENV} ${GENCODEFILTERPY} --gx ${THIS_GXFILE} \
                                    --donors ${DONORFILE} --dataset gtex \
                                    --gtf ${GENCODEFILE} \
                                    --biotype ${GXSELECTION} # --out ${TISSUEGXFILE} 
    fi
done
