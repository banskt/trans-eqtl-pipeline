#!/bin/sh
#SBATCH -p hh
#SBATCH -t 2-00:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load gcc

TMPDIR="/local/tmp/_JOB_NAME"
mkdir -p $TMPDIR && export $TMPDIR

rm -rf /tmp/*

VCFTOOLS=_VCFTOOL_
GTFILE=_GTFILE_
OUTFILE=_OUTFILE_
POP1=_POP1_
POP2=_POP2_
POP3=_POP3_

EXTRAPOP="--weir-fst-pop ${POP3}"
if [ "${POP3}" = "" ]; then
    EXTRAPOP=""
fi

${VCFTOOLS} --gzvcf ${GTFILE} \
            --out ${OUTFILE} \
            --weir-fst-pop ${POP1} \
            --weir-fst-pop ${POP2} ${EXTRAPOP}

cat ${OUTFILE}.weir.fst | sed 's/chr//g' > ${OUTFILE}.fix.weir.fst
rm -f ${OUTFILE}.weir.fst
mv ${OUTFILE}.fix.weir.fst ${OUTFILE}.weir.fst
rm -rf /tmp/*
