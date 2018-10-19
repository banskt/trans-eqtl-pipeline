#!/bin/bash

LDSTORE=$1; shift
PLOTFILE=$1; shift
CHRM=$1; shift
NT=$1; shift
BGENFILE=$1; shift
RESFILE=$1; shift
OUTDIR=$1; shift

module load plink/2.00
source activate py36

for n in $(seq 1 $#); do

    GENENAME=$1; shift

    CURDIR=`pwd`
    cd ${OUTDIR}
        cat ${GENENAME}.txt | cut -d" " -f1 > ${GENENAME}.snplist
        plink2 --bgen ${BGENFILE} --oxford-single-chr ${CHRM} --extract  ${GENENAME}.snplist --export bgen-1.2 --out ${GENENAME}
        ${LDSTORE} --bgen ${GENENAME}.bgen --bcor ${GENENAME}.bcor --n-threads ${NT}
        ${LDSTORE} --bcor ${GENENAME}.bcor --merge ${NT}
        rm -rf ${GENENAME}.bcor_*
        ${LDSTORE} --bcor ${GENENAME}.bcor --matrix ${GENENAME}.ld
        ${LDSTORE} --bcor ${GENENAME}.bcor --meta ${GENENAME}.meta
        rm -rf ${GENENAME}.bgen ${GENENAME}.bcor ${GENENAME}.snplist ${GENENAME}.log ${GENENAME}.sample
    cd ${CURDIR}
    
    python ${PLOTFILE} --gene ${GENENAME} --outdir ${OUTDIR} --meqtl ${RESFILE}

done
