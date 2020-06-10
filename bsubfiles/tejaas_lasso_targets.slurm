#!/bin/sh
#SBATCH -p hh
#SBATCH -t 0-6:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load intel/mpi/64/2017/2.174

# source $HOME/miniconda3/envs/py36/bin/activate py36

PYENV=_PY_ENV_
PYRRLASSO=_RLASSO_
SAMPLEFILE=_SAMFIL_
GXFILE=_GXFILE_
GTFILE=_VCFFIL_
GENEINFOFILE=_GENINF_
SPECIFIC_OUTDIR=_OUTDIR_
BASEDIR=_BASEDR_
CHRM=_CHRM_N_
TSHORT=_TISSUE_
K=_KNN_K_
EXPRCORRFILE=_GXCORR_

${PYENV} ${PYRRLASSO}   --tissue   ${TSHORT} \
                        --outdir   ${SPECIFIC_OUTDIR} \
                        --gtf      ${GENEINFOFILE} \
                        --basedir  ${BASEDIR} \
                        --vcf      ${GTFILE} \
                        --fam      ${SAMPLEFILE} \
                        --gxfile   ${GXFILE} \
                        --gxcorr   ${EXPRCORRFILE} \
                        --knn      ${K} \
                        --chrm     ${CHRM}
                        

### ${MPYTHON} ${RUN_PATH} --vcf ${GENOFILE} --gx ${GXPRFILE} --gtf ${GENEINFO} --method ${TJMETHOD} --null ${NULLMODL} --outprefix ${OUTPRFIX} --include-SNPs ${INCSTRNG} --psnpthres ${SNPTHRES} --pgenethres ${GENTHRES} --prior-sigma ${SBETA} --chrom ${CHROM} ${EXTRAFLAGS}
