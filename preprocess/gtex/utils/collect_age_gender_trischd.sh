#!/bin/bash

# SRCSUBJCT
# COVARDIR
# PREPROC_UTILSDIR

## collect subjectID, age, trischd from phenotypefile
grep -v -P "^#" ${SRCSUBJCT} | cut -f 2,5,15 | grep "SUBJID\|GTEX" > ${COVARDIR}/tmp_gender_age_trischd.txt

## transpose 
awk -f ${PREPROC_UTILSDIR}/transpose.awk ${COVARDIR}/tmp_gender_age_trischd.txt > ${COVARDIR}/cov_gender_age_trischd.txt
rm -rf ${COVARDIR}/tmp_gender_age_trischd.txt
