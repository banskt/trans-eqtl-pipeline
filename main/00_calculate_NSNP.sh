#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars

source ${DATALOAD}

#CHRM_NTOT_FILE="$( dirname ${GENO_FMT} )/ntot_per_chromosome.txt"

#if [ -f ${CHRM_NTOT_FILE} ]; then rm -f ${CHRM_NTOT_FILE}; touch ${CHRM_NTOT_FILE}; fi

for CHRM in ${CHRNUMS}; do
    GENOTYPEFILE=${GENO_FMT/\[CHRM\]/${CHRM}}
    SNPLISTFILE="$( dirname ${GENOTYPEFILE} )/$( basename ${GENOTYPEFILE} .vcf.gz ).snplist"
    echo $SNPLISTFILE
    zcat ${GENOTYPEFILE} | cut -f1-3 | sed '/^\s*#/d;/^\s*$/d' | cut -f3 > ${SNPLISTFILE}
    #NTOT=$( zcat ${GENOTYPEFILE} | sed '/^\s*#/d;/^\s*$/d' | wc -l )
    #echo "${CHRM}: ${NTOT} SNPs"
    #echo "${CHRM} ${NTOT}" >> ${CHRM_NTOT_FILE}
done

unset_vars ${CONFIGFILE}
unset_vars PATHS
