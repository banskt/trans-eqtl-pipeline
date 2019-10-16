#!/bin/sh

for CHROM in `seq 1 22`;
do

RUNTIME="48:00"

INPUTDIR="/cbscratch/franco/datasets"
GTFILE="${INPUTDIR}/cardiogenics/genotypes/CG_${CHROM}.imputed.gz"
DONORS="${INPUTDIR}/cardiogenics/genotypes/CG.sample"

OUTDIR="${INPUTDIR}/cardiogenics/genotypes/prefiltered_dosages"
OUTFILE="${OUTDIR}/CG_dosages_filtered_${CHROM}.imputed.gz"

DATASET="cardiogenics"  # filtered genotypes have gtex format
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
