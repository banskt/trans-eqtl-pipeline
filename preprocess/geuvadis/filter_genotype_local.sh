#!/bin/bash -e

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi
CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}

MAF_THRES="0.01"

OUTVCFDIR="${INVCFDIR}/${MAF_THRES}"
if [ ! -e $OUTVCFDIR ]; then
    mkdir -p $OUTVCFDIR
fi

INVCF="${INVCFDIR}/${VCFIBASE}.vcf.gz"
OUTVCF="${OUTVCFDIR}/${VCFOBASE}.vcf.gz"

for CHRM in `seq 1 21`; do
    CHRM_INVCF=${INVCF/\[CHRM\]/${CHRM}}
    CHRM_OUTVCF=${OUTVCF/\[CHRM\]/${CHRM}}
    CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

    CHRM_OUTVCF=${CHRM_OUTVCF/\[MAFTHRES\]/${MAF_THRES}}

    $PYENV $PYFILTER --in $CHRM_INVCF --out $CHRM_OUTVCF --maf $MAF_THRES &
done