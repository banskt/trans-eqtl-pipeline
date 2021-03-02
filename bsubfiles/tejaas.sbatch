#!/bin/bash
#SBATCH -A cramer
#SBATCH -p em
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME_%a.out
#SBATCH -e _JOB_NAME_%a.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
# module load openmpi/intel/64/1.10.7

export OMP_NUM_THREADS=1
MPYTHON=_PYT_ENV_

RUN_PATH=_TJS_BINR
GENOFILE=_GT_FILE_
GXPRFILE=_EXPR_FL_
GXCORRFL=_GX_CRFL_
SAMPLEFL=_SAMP_FL_
GENEINFO=_GEN_POSF

TJMETHOD=_TJ_METHD
NULLMODL=_NULL_MDL
INDEX=$( echo ${SLURM_ARRAY_TASK_ID} | awk '{printf "%03d", $1}' )
OUTPRFIX=_OUT_PRFX${INDEX}
SNPTHRES=_SNP_CUT_
GENTHRES=_GEN_CUT_
SBETA=_SIG_BETA
CHROM=_CHRM_NUM
MAFFILE=_MAF_FILE
EXTRAFLAGS=_EXT_FLAG

if [ "${NULLMODL}" = "maf" ]; then
    EXTRAFLAGS="$EXTRAFLAGS --maf-file ${MAFFILE}"
fi


NMAX=_NUM_MAX_
NTOT=_NUM_TOT_
STARTSNP=$(( NMAX * SLURM_ARRAY_TASK_ID + 1 ))
ENDSNP=$(( NMAX * (SLURM_ARRAY_TASK_ID + 1) ))
if [ $ENDSNP -gt $NTOT ]; then
    ENDSNP=${NTOT}
fi
INCSTRNG="${STARTSNP}:${ENDSNP}"

if [ "${NULLMODL}" = "maf" ]; then
    EXTRAFLAGS="$EXTRAFLAGS --maf-file ${MAFFILE}"
fi



mpirun -n 8 ${MPYTHON} ${RUN_PATH} --vcf          ${GENOFILE} \
                                   --gx           ${GXPRFILE} \
                                   --gxcorr       ${GXCORRFL} \
                                   --fam          ${SAMPLEFL} \
                                   --gtf          ${GENEINFO} \
                                   --method       ${TJMETHOD} \
                                   --null         ${NULLMODL} \
                                   --outprefix    ${OUTPRFIX} \
                                   --include-SNPs ${INCSTRNG} \
                                   --psnpthres    ${SNPTHRES} \
                                   --pgenethres   ${GENTHRES} \
                                   --prior-sigma  ${SBETA} \
                                   --chrom        ${CHROM} \
                                   ${EXTRAFLAGS}

### ${MPYTHON} ${RUN_PATH} --vcf ${GENOFILE} --gx ${GXPRFILE} --gxcorr ${GXCORRFL} --gtf ${GENEINFO} --method ${TJMETHOD} --null ${NULLMODL} --outprefix ${OUTPRFIX} --include-SNPs ${INCSTRNG} --psnpthres ${SNPTHRES} --pgenethres ${GENTHRES} --prior-sigma ${SBETA} --chrom ${CHROM} ${EXTRAFLAGS}
