#!/bin/sh -e

# PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
# PYTPMS="${PWD}/scripts/calculated_TPMs.py"

${PYENV} ${PYTPMS} --in ${SRCREAD} --out ${OUTTPM} --gtf ${GTFFILE}