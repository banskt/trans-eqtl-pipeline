#!/bin/bash
GTEXCOV="${SRCCOVARS}/${TBASE}_Analysis.v6p.covariates.txt"
echo "Processing Covariates for Tissue: $TFULL"
# Selects only genotype PCs, sex and platform
grep -v -i "inferred" ${GTEXCOV} > ${COVARS}

sleep 5

${PYENV} ${COMPILEAGECOVPY} --input ${COVARS} --age ${AGE_COVARIATE_FILE} --output ${COVARS_AGE}
