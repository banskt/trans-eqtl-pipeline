#!/bin/bash
# mkdir -p ${RPKMOUTDIR};
echo "Processing Tissue: $TFULL"
mkdir -p $RPKMOUTDIR
RPKMFILE="${RPKMOUTDIR}/${TSHORT}_rpkm.gct"

# Select tissue-specific rpkms
OUTPREFIXFILE="${RPKMOUTDIR}/${TSHORT}"
${PYENV} ${SELECTSAMPLEPY} --rpkm $SRCRPKM --counts $SRCREAD --output $OUTPREFIXFILE --tissue="$TFULL" --pheno $SRCPHENO

# Normalize
RPKMFILE="${RPKMOUTDIR}/${TSHORT}_rpkm.gct"
COUNTSFILE="${RPKMOUTDIR}/${TSHORT}_counts.gct"
${PYENV} ${GTEXNORMALIZEPY} --rpkm $RPKMFILE --counts $COUNTSFILE --tissue $TSHORT --donors $DONORFILE --outdir $GXOUTDIR


