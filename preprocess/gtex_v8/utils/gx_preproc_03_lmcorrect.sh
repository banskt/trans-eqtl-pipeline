#!/bin/bash -e

TS="${TSHORT}"
if [ '${bUsePub}' == "true" ]; then
    TS="${TSHORT}_pub"
fi

${PYENV} ${CORRECTPY} --input-dir ${GXOUTDIR} --tissue ${TS} --cov ${COVARS_AGE} --outdir ${GXOUTDIR}
