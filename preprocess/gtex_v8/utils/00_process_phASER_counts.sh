#!/bin/sh -e

PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
PYPHASER="${PWD}/scripts/process_phASER_counts.py"

INMTX="/cbscratch/franco/datasets/gtex_v8/expression/phASER_GTEx_v8_matrix.txt.gz"
OUTMTX="/cbscratch/franco/datasets/gtex_v8/expression/phASER_GTEx_v8_matrix_collapsed_counts_new.txt.gz"

${PYENV} ${PYPHASER} --in ${INMTX} --out ${OUTMTX}