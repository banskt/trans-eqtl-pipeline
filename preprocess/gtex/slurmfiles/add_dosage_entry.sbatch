#!/bin/bash
#SBATCH -p medium
#SBATCH -t 2-0:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -C scratch
#SBATCH -J add_DSentry
#SBATCH -o add_DSentry_%A_%a.out
#SBATCH -e add_DSentry_%A_%a.err

#CHROM=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
CHROM=(1 2)
DS2VCF="/usr/users/sbanerj/trans-eQTL/dev-pipeline/preprocess/external/DS2VCF/bin/add_dosage"
CHRNUM=${CHROM[$SLURM_ARRAY_TASK_ID]}
INFILE="/scratch/sbanerj/trans-eqtl/input/gtex_v8/genotype/all_samples/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_chr${CHRNUM}.vcf.gz"
OUTFILE="/usr/users/sbanerj/gtex_v8/genotype/all_samples/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_withDS_chr${CHRNUM}.vcf.gz"
${DS2VCF} ${INFILE} ${OUTFILE}
