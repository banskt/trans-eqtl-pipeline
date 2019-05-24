#!/bin/bash
mkdir -p ${RPKMOUTDIR};
mkdir -p ${NORMOUTDIR}
echo "Processing Tissue: $TFULL"
RPKMFILE="${RPKMOUTDIR}/${TSHORT}_rpkm.gct"
${PYENV} ${SELECTSAMPLEPY} --input ${SRCRPKM} --output ${RPKMFILE} --tissue="$TFULL" --pheno ${SRCPHENO}
${PYENV} ${GTEXNORMALIZEPY} --rpkm ${RPKMFILE} --counts ${SRCREAD} --tissue ${TSHORT} --donors ${DONORFILE} --outdir ${NORMOUTDIR}
