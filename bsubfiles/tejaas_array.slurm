#!/bin/sh
#SBATCH -p hh
#SBATCH -t 0-12:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J _JOB_NAME_%a
#SBATCH -o _JOB_NAME_%a.out
#SBATCH -e _JOB_NAME_%a.err

#module load intel/compiler/64/2020/19.1.2
#module load intel/mkl/64/2020/2.254
#module load intel/mpi/64/2020/2.254

module load intel/2020.4 intel-parallel-studio/cluster.2020.4 # intel-mpi/2019.9.304
module load openmpi/4.0.5


# source $HOME/miniconda3/envs/py36/bin/activate py36

RUN_PATH=_TJS_BINR
GENOFILE=_GT_FILE_
SAMPFILE=_SAM_FILE
GXPRFILE=_EXPR_FL_
GENEINFO=_GEN_POSF

TJMETHOD=_TJ_METHD
NULLMODL=_NULL_MDL
INDEX=`echo ${SLURM_ARRAY_TASK_ID} | awk '{printf "%03d", $1}'`
OUTPRFIX=_OUT_PRFX${INDEX}
SNPTHRES=_SNP_CUT_
GENTHRES=_GEN_CUT_
SBETA=_SIG_BETA
CHROM=_CHRM_NUM
MAFFILE=_MAF_FILE
EXTRAFLAGS=_EXT_FLAG

NMAX=_NUM_MAX_
NTOT=_NUM_TOT_
startsnp=$(( NMAX * SLURM_ARRAY_TASK_ID + 1 ))
endsnp=$(( NMAX * (SLURM_ARRAY_TASK_ID + 1) ))
if [ $endsnp -gt $NTOT ]; then
    endsnp=${NTOT}
fi
INCSTRNG="${startsnp}:${endsnp}"

if [ "${NULLMODL}" = "maf" ]; then
    EXTRAFLAGS="$EXTRAFLAGS --maf-file ${MAFFILE}"
fi

# mpirun -n 8 ${RUN_PATH} --oxf          ${GENOFILE} \
#                         --fam          ${SAMPFILE} \
#                         --gx           ${GXPRFILE} \
#                         --gtf          ${GENEINFO} \
#                         --method       ${TJMETHOD} \
#                         --null         ${NULLMODL} \
#                         --outprefix    ${OUTPRFIX} \
#                         --include-SNPs ${INCSTRNG} \
#                         --psnpthres    ${SNPTHRES} \
#                         --pgenethres   ${GENTHRES} \
#                         --prior-sigma  ${SBETA} \
#                         --chrom        ${CHROM} \
#                         ${EXTRAFLAGS}

mpirun -n 8 ${RUN_PATH} --vcf          ${GENOFILE} \
                        --fam          ${SAMPFILE} \
                        --gx           ${GXPRFILE} \
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

### ${MPYTHON} ${RUN_PATH} --vcf ${GENOFILE} --gx ${GXPRFILE} --gtf ${GENEINFO} --method ${TJMETHOD} --null ${NULLMODL} --outprefix ${OUTPRFIX} --include-SNPs ${INCSTRNG} --psnpthres ${SNPTHRES} --pgenethres ${GENTHRES} --prior-sigma ${SBETA} --chrom ${CHROM} ${EXTRAFLAGS}
