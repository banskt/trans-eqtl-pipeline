#!/bin/sh
#SBATCH -p hh
#SBATCH -t 0-12:00:00
#SBATCH -n _NCORE_
#SBATCH -N 1
#SBATCH -J compreq_chr1_snp_NSNPS__gene_NGENE__samp_NSAMP__ncore_NCORE_
#SBATCH -o logs/compreq_chr1_snp_NSNPS__gene_NGENE__samp_NSAMP__ncore_NCORE_.out
#SBATCH -e logs/compreq_chr1_snp_NSNPS__gene_NGENE__samp_NSAMP__ncore_NCORE_.err
#SBATCH -x hh036,hh037,hh038,hh039,hh040


## G Genes: 15673 as_tpms_qcfilter.txt.protein_coding_lncRNA_filtered


module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load intel/mpi/64/2017/2.174

RUN_PATH=/usr/users/fsimone/tejaas/bin/tejaas
GENOFILE=/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs_SHAPEIT2/0.01/GTEX_v8_2020-02-21_WGS_838Indiv_Freeze.SHAPEIT2_phased_NoMissingGT_SNPfilter_MAF0.01_chr1.vcf.gz
SAMPFILE=/usr/users/fsimone/trans-eqtl-pipeline/comp_req/samples/gtex_v8_as_N_NSAMP_.txt
#GXPRFILE=/cbscratch/franco/trans-eqtl/new_preprocess_feb2020_freeze/gtex_v8/expression/tpms/as_tpms_qcfilter.txt.protein_coding_lncRNA_filtered
GXPRFILE=/usr/users/fsimone/trans-eqtl-pipeline/comp_req/expr/as_tpms_qcfilter_N_NGENE_.txt
GENEINFO=/cbscratch/franco/datasets/GENCODE/gencode.v26.annotation.gtf.gz

TJMETHOD=jpa-rr
NULLMODL=perm
INDEX=`echo ${SLURM_ARRAY_TASK_ID} | awk '{printf "%03d", $1}'`
OUTPRFIX=/cbscratch/franco/trans-eqtl/dev-pipeline/tejaas_comp_req/raw/gtex_v8-as/tejaas/permnull_sb0.1_knn30/chr1/chunk${INDEX}
SNPTHRES=0.0000005
GENTHRES=0.01
SBETA=0.1
CHROM=1
MAFFILE=None
EXTRAFLAGS="--dosage --gxcorr /usr/users/fsimone/trans-eqtl-pipeline/comp_req/expr/as_tpms_cclm_N_NGENE_.txt --knn 30 --cismask"

NMAX=_NSNPS_
NTOT=629624
startsnp=0     #$(( NMAX * SLURM_ARRAY_TASK_ID + 1 ))
endsnp=$NMAX #$(( NMAX * (SLURM_ARRAY_TASK_ID + 1) ))
if [ $endsnp -gt $NTOT ]; then
    endsnp=${NTOT}
fi
INCSTRNG="${startsnp}:${endsnp}"

mpirun -n _NCORE_ ${RUN_PATH} --vcf          ${GENOFILE} \
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
