#!/bin/sh
#SBATCH -p hh
#SBATCH -t 0-12:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

# module load R/3.4.1

SCRIPT=_MEQTL_R_
GTFILE=_GT_FILE_
SAMPLEFILE=_SAM_FILE
GENPOSFILE=_GEN_POSF
EXPRFILE=_EXPR_FL_
DATATYPE=_MD_TYPE_
CISPVAL=_PTHR_CIS
TRANSPVAL=_PTHR_TRA
OUTFILECIS=_OUTF_CIS
OUTFILETRANS=_OUTF_TRA
MATRIXEQTLMODEL=_MEQTLMOD
EXTRAFLAGS=_EXT_FLAG


Rscript ${SCRIPT} --genotype ${GTFILE} \
                  --donors ${SAMPLEFILE} \
                  --geneinfo ${GENPOSFILE} \
                  --datatype ${DATATYPE} \
                  --expression ${EXPRFILE} \
                  --pvalcis ${CISPVAL} \
                  --pvaltrans ${TRANSPVAL} \
                  --outfilecis ${OUTFILECIS} \
                  --outfiletrans ${OUTFILETRANS} \
                  --model ${MATRIXEQTLMODEL} \
                  ${EXTRAFLAGS}

PYENV=_PYENV_
PYMQTL=_PYMQTL_
POSFILE=_POSFILE_

${PYENV} ${PYMQTL} --input ${OUTFILECIS} --posfile ${POSFILE}
${PYENV} ${PYMQTL} --input ${OUTFILETRANS} --posfile ${POSFILE}