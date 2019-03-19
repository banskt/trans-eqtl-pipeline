#!/bin/bash -e

mkdir -p $NORMOUTDIR

echo "Processing Tisuee: $TFULL"
# Select tissue-specific rpkms
RPKMFILE="${RPKMOUTDIR}/${TSHORT}_rpkm.gct"
$PYENV ${SCRIPTDIR}/select_samples_from_tissue.py --input $SRCRPKM --output $RPKMFILE --tissue="$TFULL" --pheno $SRCPHENO

# Normalize
$PYENV ${SCRIPTDIR}/do_expression_normalization.py --rpkm $RPKMFILE --counts $SRCREAD --tissue $TSHORT --donors $DONORFILE --outdir $NORMOUTDIR