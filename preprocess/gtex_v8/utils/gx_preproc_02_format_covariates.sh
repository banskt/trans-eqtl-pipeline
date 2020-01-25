#!/bin/bash -e
GTEXCOV="${SRCCOVARS}/${TBASE}.v8.covariates.txt"
echo "Processing Covariates for Tissue: $TFULL"
# Selects only genotype PCs, sex and platform

## this also deletes PC1-5 of the genotype
grep -v -i "inferred" ${GTEXCOV} |grep -v -i -P "pc\d" > ${COVARS_NOPC}
grep -v -i "inferred" ${GTEXCOV} > ${COVARS}

sleep 5

${PYENV} ${COMPILEAGECOVPY} --input ${COVARS} --age ${AGE_COVARIATE_FILE} --output ${COVARS_AGE}
${PYENV} ${COMPILEAGECOVPY} --input ${COVARS_NOPC} --age ${AGE_COVARIATE_FILE} --output ${COVARS_AGE_NOPC}
