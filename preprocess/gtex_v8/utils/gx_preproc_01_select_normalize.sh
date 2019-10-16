#!/bin/bash -e
echo "Processing Tissue: $TFULL"
mkdir -p $TPMOUTDIR

TS="${TSHORT}"
if [ '${bUsePub}' == "true" ]; then
    TS="${TSHORT}_pub"
fi

# Select tissue-specific tpms
OUTPREFIXFILE="${TPMOUTDIR}/${TS}"
echo ${PYENV} ${SELECTSAMPLEPY} --rpkm $SRCTPM --counts $SRCREAD --output $OUTPREFIXFILE --tissue="$TFULL" --pheno $SRCPHENO
${PYENV} ${SELECTSAMPLEPY} --rpkm $SRCTPM --counts $SRCREAD --output $OUTPREFIXFILE --tissue="$TFULL" --pheno $SRCPHENO

# Normalize
TPMFILE="${TPMOUTDIR}/${TS}_tpm.gct"
COUNTSFILE="${TPMOUTDIR}/${TS}_counts.gct"
echo ${PYENV} ${GTEXNORMALIZEPY} --rpkm $TPMFILE --counts $COUNTSFILE --tissue $TS --donors $DONORFILE --outdir $GXOUTDIR
${PYENV} ${GTEXNORMALIZEPY} --rpkm $TPMFILE --counts $COUNTSFILE --tissue $TS --donors $DONORFILE --outdir $GXOUTDIR
