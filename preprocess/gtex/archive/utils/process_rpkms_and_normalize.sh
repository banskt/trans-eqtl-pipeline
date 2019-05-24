#!/bin/bash -e

mkdir -p $NORMOUTDIR

echo "Processing Tissue: $TFULL"
# Select tissue-specific rpkms
OUTPREFIXFILE="${RPKMOUTDIR}/${TSHORT}"
$PYENV ${SCRIPTDIR}/select_samples_from_tissue.py --rpkm $SRCRPKM --counts $SRCREAD --output $OUTPREFIXFILE --tissue="$TFULL" --pheno $SRCPHENO

# Normalize
RPKMFILE="${RPKMOUTDIR}/${TSHORT}_rpkm.gct"
COUNTSFILE="${RPKMOUTDIR}/${TSHORT}_counts.gct"
$PYENV ${SCRIPTDIR}/do_expression_normalization.py --rpkm $RPKMFILE --counts $COUNTSFILE --tissue $TSHORT --donors $DONORFILE --outdir $NORMOUTDIR