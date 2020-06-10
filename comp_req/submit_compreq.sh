#!/bin/bash

NGENEMAX=15673 #(+1, for the header)
NGENES="5000 10000 ${NGENEMAX}"

MAXNSNP=629624
NSNPS="100 1000 5000 10000 20000 40000 60000 80000 100000 120000"

NSAMPMAX="581" #as samples with genotype
NSAMPLES="200 300 400 500 ${NSAMPMAX}"

for NSAMP in $NSAMPLES; do
      for NGENE in $NGENES; do
            for NSNP in $NSNPS; do

                  slurmfile="compreq_snp${NSNP}_gene${NGENE}_samp${NSAMP}.slurm"
                  sed "s|_NSNPS_|${NSNP}|g;
                       s|_NGENE_|${NGENE}|g;
                       s|_NSAMP_|${NSAMP}|g;
                       " compreq_basejob.slurm > "jobsubs/${slurmfile}"

                  sbatch "jobsubs/${slurmfile}"
            done;
      done;
done;
