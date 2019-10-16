#!/bin/bash -e
GTEXCOV="${SRCCOVARS}/${TBASE}.v8.covariates.txt"
echo "Processing Covariates for Tissue: $TFULL"
# Selects only genotype PCs, sex and platform
grep -v -i "inferred" ${GTEXCOV} > ${COVARS}

sleep 5

${PYENV} ${COMPILEAGECOVPY} --input ${COVARS} --age ${AGE_COVARIATE_FILE} --output ${COVARS_AGE}
