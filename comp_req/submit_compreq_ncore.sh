#!/bin/bash

NGENEMAX=15673 #(+1, for the header)
NGENES="${NGENEMAX}"

MAXNSNP=629624
NSNPS="10000 20000 40000"

NSAMPMAX="581" #as samples with genotype
NSAMPLES="${NSAMPMAX}"

NCORES="1 2 4 8 16"

for NCORE in $NCORES; do
    for NSAMP in $NSAMPLES; do
        for NGENE in $NGENES; do
                for NSNP in $NSNPS; do

                    slurmfile="compreq_snp${NSNP}_gene${NGENE}_samp${NSAMP}_ncore${NCORE}.slurm"
                    sed "s|_NSNPS_|${NSNP}|g;
                        s|_NGENE_|${NGENE}|g;
                        s|_NSAMP_|${NSAMP}|g;
                        s|_NCORE_|${NCORE}|g;
                        " compreq_basejob_ncore.slurm > "jobsubs/${slurmfile}"

                    sbatch "jobsubs/${slurmfile}"
                done;
        done;
    done;
done;
