#!/bin/bash
mkdir -p ${PEEROUTDIR};

EXPR_TYPES="qn tmm tpms"

## PEER correction on normalized gene expression (without using any covariates)
for EXPR_TYPE in ${EXPR_TYPES}; do
    for NPEER in ${GXNPEERS}; do
        THIS_PEEROUTDIR="${PEEROUTDIR}/${EXPR_TYPE}"
        INFILE="${EXPR_TYPE}/${TSHORT}_${EXPR_TYPE}.txt"
        if [ ! -d "${THIS_PEEROUTDIR}" ]; then mkdir -p ${THIS_PEEROUTDIR}; fi
        Rscript ${PEERSCRIPT_R} ${NORMOUTFILE} "${TSHORT}_${NPEER}" --n ${NPEER} --output_dir ${THIS_PEEROUTDIR}
    done
done

for NPEER in ${GXNPEERS}; do
    THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
    if [ ! -d "${THIS_PEEROUTDIR}" ]; then mkdir -p ${THIS_PEEROUTDIR} ; fi
    Rscript ${PEERSCRIPT_R} ${NORMOUTFILE} "${TSHORT}_${NPEER}" --covar ${COVARS_AGE} --n ${NPEER} --output_dir ${THIS_PEEROUTDIR}
done    


# ## PEER correction on normalized gene expression (using other covariates)
# for bAgeCorrect in ${GXAGECORR}; do
#     for NPEER in ${GXNPEERS}; do
#         if [ "${bAgeCorrect}" = "true" ]; then
#             THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
#             if [ ! -d "${THIS_PEEROUTDIR}" ]; then mkdir -p ${THIS_PEEROUTDIR} ; fi
#             Rscript ${PEERSCRIPT_R} ${NORMOUTFILE} "${TSHORT}_${NPEER}" --covar ${COVARS_AGE} --n ${NPEER} --output_dir ${THIS_PEEROUTDIR}
#         else
#             THIS_PEEROUTDIR="${PEEROUTDIR}/covar"
#             if [ ! -d "${THIS_PEEROUTDIR}" ]; then mkdir -p ${THIS_PEEROUTDIR}; fi
#             Rscript ${PEERSCRIPT_R} ${NORMOUTFILE} "${TSHORT}_${NPEER}" --covar ${COVARS} --n ${NPEER} --output_dir ${THIS_PEEROUTDIR}
#         fi
#     done
# done
