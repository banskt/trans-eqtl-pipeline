#!/bin/sh

INPUTDIR="/cbscratch/franco/datasets/gtex/genotypes"
OUTDIR="${INPUTDIR}/prefiltered_dosages"
mkdir -p $OUTDIR

for CHROM in `seq 1 22`;
do

RUNTIME="48:00"

GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
DONORS="${INPUTDIR}/donor_ids.fam"


OUTFILE="${OUTDIR}/GTEx_450Indiv_filtered_chr${CHROM}.gz"

DATASET="gtex"  # filtered genotypes have gtex format
JOBNAME="filter_chr${CHROM}"
PYTHON="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"

mkdir -p $OUTDIR

echo "Filtering CHR${CHROM} "
bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch\
          -J ${JOBNAME} \
          -o ${OUTDIR}/${JOBNAME}.log \
          -e ${OUTDIR}/${JOBNAME}.err \
          $PYTHON filter_chromosomes.py --oxf ${GTFILE} \
                                        --chrom ${CHROM}  \
                                        --fam ${DONORS} \
                                        --out ${OUTFILE} \
                                        --dset ${DATASET}

done;
