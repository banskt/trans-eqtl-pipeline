#!/bin/bash
# mkdir -p ${LMOUTDIR};
# for bAgeCorrect in ${GXAGECORR}; do
#     if [ "${bAgeCorrect}" = "true" ]; then
#         if [ -e ${NORMOUTFILE} ] && [ -e ${COVARS} ] && [ -e ${COVARS_AGE} ]; then
#             Rscript ${LMCORR_R} ${NORMOUTFILE} ${LMOUTFILE_AGE} ${COVARS_AGE}
#         else
#             echo "${TFULL} covariates or ${NORMOUTFILE} not found"
#         fi
#     else
#         if [ -e ${NORMOUTFILE} ] && [ -e ${COVARS} ]; then
#             Rscript ${LMCORR_R} ${NORMOUTFILE} ${LMOUTFILE} ${COVARS}
#         else
#             echo "${TFULL} covariates or ${NORMOUTFILE} not found"
#         fi
#     fi
# done

${PYENV} ${CORRECTPY} --input-dir ${GXOUTDIR} --tissue ${TSHORT} --cov ${COVARS_AGE} --outdir ${GXOUTDIR}
