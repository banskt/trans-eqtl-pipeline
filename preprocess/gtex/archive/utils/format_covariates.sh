#!/bin/bash -e

# Obtain Covariates from GTEX
GTEXCOV="${SRCCOVARS}/${TBASE}_Analysis.v6p.covariates.txt"
COVARS="${COVOUTDIR}/${TSHORT}_nopeer_covariates.txt"
COVARS_AGE="${COVOUTDIR}/${TSHORT}_nopeer_covariates_w_age.txt"

echo "Processing Covariates for Tissue: $TFULL"
# Selects only genotype PCs, sex and platform
grep -v -i "inferred" $GTEXCOV > $COVARS
$PYENV ${SCRIPTDIR}/compile_age_covariate.py --input $COVARS --age $AGE_COVARIATE --output $COVARS_AGE
